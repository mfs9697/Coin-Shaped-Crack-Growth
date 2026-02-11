function run_Notch3D02()
%RUN_NOTCH3D02  Step2-style driver for Notch3D case (no globals).
%
% Implements:
%   - kdt==0: elastic initialization via F0c_core
%   - kdt>0: dt selection from scalar equation Sigma(dt)=sig_target
%   - state update: U, S-history, Fpsig, time vector
%
% Requirements:
%   cfg_Notch3D02.m must define at least:
%       C.nt, C.p, C.sig_target, C.dt_min, C.dt_max, C.dt_guess,
%       C.use_fzero, C.bracketGrow,
%       C.solve.fsolve_F0c (optimoptions for fsolve),
%       C.opt_dt_fsolve (optional if not using fzero),
%       C.opt_fzero (optimset) (optional)
%   init_Notch3D02.m must return ctx,state with the fields listed at bottom.

clc;

%% Paths
repoRoot = fileparts(fileparts(mfilename('fullpath'))); % .../cases
projRoot = fileparts(repoRoot);                         % project root
addpath(genpath(projRoot));

%% Config + init
C = cfg_Notch3D02();
[ctx, state] = init_Notch3D02(C);

fprintf('ndof=%d, nelem=%d, nnodes=%d, nip3=%d, M=%d\n', ...
    ctx.ndof, ctx.nelem, ctx.nnodes, ctx.nip3, ctx.M);
fprintf('Initial ci=%d\n\n', state.ci);

%% Main loop
for kdt = 0:C.nt
    state.kdt = kdt;

    if kdt == 0
        % Elastic init (matches legacy: K(Eo), then fsolve(F0c), then kD)
        [ctx, state] = do_initial_elastic_step(ctx, state, C);
        % After kdt==0, update stepping rule once (legacy did this after plotting)
        state = update_control_indices(ctx, state, C);
        continue;
    end

    % ---- dt selection (Sigma(dt)=sig_target) ----
    sig_target = C.sig_target;

    dt_min = C.dt_min;
    dt_max = C.dt_max;

    sig_min = Sigma_of_dt(dt_min, ctx, state, C);
    fprintf('kdt=%d: Sigma(dt_min)=%.6e, target=%.6e\n', kdt, sig_min, sig_target);

    % If Sigma(dt_min) <= target, no admissible quasi-static dt (under this control)
    if sig_min <= sig_target
        fprintf('Stop: Sigma(dt_min) <= target (no admissible quasi-static dt)\n');
        break;
    end

    if isfield(C,'use_fzero') && C.use_fzero
        g = @(dt) Sigma_of_dt(dt, ctx, state, C) - sig_target;
        [dt, info] = dt_bracket_and_fzero(g, dt_min, dt_max, C.bracketGrow, C);
        fprintf('  bracket: dtL=%.3e gL=%.3e | dtR=%.3e gR=%.3e | %s\n', ...
            info.dtL, info.gL, info.dtR, info.gR, info.status);
        if ~isfinite(dt)
            fprintf('Stop: no root bracketed in [dt_min, dt_max].\n');
            break;
        end
    else
        % legacy-ish scalar fsolve
        if ~isfield(C,'opt_dt_fsolve') || isempty(C.opt_dt_fsolve)
            C.opt_dt_fsolve = optimoptions('fsolve','Display','off','MaxIterations',50);
        end
        dt = fsolve(@(dt) Sigma_of_dt(dt, ctx, state, C) - sig_target, C.dt_guess, C.opt_dt_fsolve);
    end

    state.dt = dt;

    % ---- accept step: compute dU + update history consistently ----
    [ctx, state, diag] = advance_one_increment(ctx, state, C); %#ok<NASGU>

    % ---- plots (optional) ----
    if isfield(C,'doPlot') && C.doPlot
        do_plots(ctx, state, C);
    end

    % ---- update kD / ci ----
    state = update_control_indices(ctx, state, C);

    if isfield(C,'pauseSec') && C.pauseSec > 0
        pause(C.pauseSec);
    end
end

fprintf('\nrun_Notch3D02 finished: %s\n', datestr(now));
end

%% ========================================================================
% kdt==0 elastic step (legacy-faithful via F0c_core)
% ========================================================================
function [ctx, state] = do_initial_elastic_step(ctx, state, C)

% Elastic stiffness
if ~isfield(ctx,'Eo')
    error('do_initial_elastic_step: ctx.Eo missing (elastic modulus matrix).');
end
ctx.K = StifFnd(ctx.Eo, ctx);

% Apply Dirichlet constraints (row/col pinning)
fix = ctx.fixvars(:);
ctx.K(fix,:) = 0; ctx.K(:,fix) = 0;
ctx.K(sub2ind(size(ctx.K), fix, fix)) = 1;

% Solve du under fixed external traction level = state.sig (set by init)
state.Fpsig = zeros(ctx.ndof,1);
state.Ftsig = zeros(ctx.ndof,1);

fun = @(du) F0c_core(du, ctx, state);

if ~isfield(C,'solve') || ~isfield(C.solve,'fsolve_F0c')
    error('do_initial_elastic_step: C.solve.fsolve_F0c missing.');
end

du = fsolve(fun, state.dU0, C.solve.fsolve_F0c);
du = full(du);

state.dU = du;
state.dU0 = du;              % good initial guess for later
state.U  = state.U + du;     % accept the elastic increment

% Initial kD from first control DOF
state.kD = 2*du(state.vc(1)) / ctx.Dym;

% Initialize time
state.t1(1) = 0;

fprintf('kdt=0: sig_ext=%.6e | kD=%.6e\n', state.sig, state.kD);
end

%% ========================================================================
% Sigma(dt) probe (non-mutating for caller)
% ========================================================================
function sig = Sigma_of_dt(dt, ctx, state, C)
    [sig, ~] = Ucurr_core(dt, C, ctx, state);
end

%% ========================================================================
% Accept one increment: call Ucurr_core then update U,S,Fpsig,time
% ========================================================================
function [ctx, state, diag] = advance_one_increment(ctx, state, C)
diag = struct();

% --- sanity ---
if ~isfield(ctx,'Nextr') || isempty(ctx.Nextr)
    error('advance_one_increment:MissingNextr', ...
        'ctx.Nextr is missing. Call [xip3,w3,Nextr]=IntegrationParameters(nip3,nelnodes,@N_3DT) in init.');
end
if ~isfield(state,'sigzip') || isempty(state.sigzip)
    state.sigzip = zeros(ctx.nelem*ctx.nip3,1);
end
if ~isfield(state,'sigz') || isempty(state.sigz)
    state.sigz = zeros(ctx.nnodes,1);
end

% 1) Get dU, K(dt), Ftsig(dt), zeta/eta(dt), and sigma_out
[sig_out, state_new] = Ucurr_core(state.dt, C, ctx, state);

% Adopt fields computed by Ucurr_core
state.zetamn = state_new.zetamn;
state.etamn  = state_new.etamn;
state.Ftsig  = state_new.Ftsig;
state.K      = state_new.K;
state.tE     = state_new.tE;
state.dU     = state_new.dU;

state.sig = sig_out;
dU = state.dU;

% 2) Update displacement
state.U   = state.U + dU;
state.dU0 = dU;  % next initial guess

% 3) Update internal history S + update sigzip (LEGACY-ORDER) + build nodal sigz
M     = ctx.M;
indU  = zeros(ctx.eldf,1);

sigz_acc = zeros(ctx.nnodes,1);
w_acc    = zeros(ctx.nnodes,1);

for el = 1:ctx.nelem
    nod = ctx.connect(:,el);
    for j = 1:ctx.nelnodes
        indU((3*j-2):(3*j)) = (3*nod(j)-2):(3*nod(j));
    end

    % loop IPs
    for ip = 1:ctx.nip3
        B3    = B(ctx.xip3(:,ip), el, ctx);
        deps  = B3 * dU(indU);                    % 6x1
        lmnip = (el-1)*ctx.nip3 + ip;

        if state.kdt > 0
            % tilda sigma from visco history (uses OLD S before update, like legacy)
            tsig = zeros(6,1);
            for m = 1:M
                tsig = tsig + sum((1 - state.zetamn(:,:,m)).*state.S(:,:,lmnip,m), 2);
            end

            % update sigma_zz at IP (store only component 3, like legacy sigzip)
            state.sigzip(lmnip) = state.sigzip(lmnip) - tsig(3) + state.tE(3,:)*deps;
        else
            % initial elastic step
            state.sigzip(lmnip) = ctx.Eo(3,:)*deps;
        end

        % now update visco history S (same as legacy)
        deps1 = repmat(deps', 6, 1);
        for m = 1:M
            state.S(:,:,lmnip,m) = state.zetamn(:,:,m).*state.S(:,:,lmnip,m) + ...
                state.etamn(:,:,m).*(1 - state.zetamn(:,:,m)).*deps1;
        end
    end

    % Extrapolate IP sigma_zz -> nodal sigma_zz for this element, then volume-average
    elIPs = (el-1)*ctx.nip3 + (1:ctx.nip3);
    sigzv = ctx.Nextr * state.sigzip(elIPs);

    w = ctx.velem(el);
    for a = 1:ctx.nelnodes
        sigz_acc(nod(a)) = sigz_acc(nod(a)) + w*sigzv(a);
        w_acc(nod(a))    = w_acc(nod(a)) + w;
    end
end

% Final nodal sigma_z: volume-weight average, scale by sc1, cap at 1
sigz = sigz_acc ./ max(w_acc, eps);
sigz = sigz ./ ctx.sc1;
sigz(sigz > 1) = 1;

state.sigz = sigz;

% 4) Update force history (legacy: Fpsig = K*dU + Fpsig - Ftsig)
state.Fpsig = state.K*dU + state.Fpsig - state.Ftsig;

% 5) Update time vector
kdt = state.kdt;
if kdt == 1
    state.t1(2) = state.dt;
elseif kdt > 1
    state.t1(kdt+1) = state.t1(kdt) + state.dt;
end

fprintf('kdt=%d: dt=%.4e | sig=%.6e | ci=%d | kD=%.6e\n', ...
    kdt, state.dt, state.sig, state.ci, state.kD);

diag.sig = state.sig;
diag.dt  = state.dt;
end


%% ========================================================================
% Plot hook (optional)
% ========================================================================
function do_plots(ctx, state, C) %#ok<INUSD>
    fact=100; 
    Mesh.Points   = ctx.coords+fact*reshape(state.U,3,ctx.nnodes);
    Mesh.Elements = ctx.connect;
    figure(1); clf(1)
    
    plot3D(Mesh,'ColorMapData',state.sigz); hold on
    plot3D(Mesh,'FaceAlpha',0,'EdgeColor',.6*[1,1,1]); hold on
    
    col=(linspace(-.1,1,12))';
    cmap=jet(length(col)-1);
    clim([col(1) col(end)]); colormap(cmap);
    colorbar('Ticks',col); hold on  
    
    view(60,-30); axis equal
end

%% ========================================================================
% Update kD/ci stepping rule (legacy-faithful)
% ========================================================================
function state = update_control_indices(ctx, state, C) %#ok<INUSL>
if ~isfield(C,'epsTol') || isempty(C.epsTol)
    C.epsTol = 1e-12;
end
e = C.epsTol;
p = C.p;

if state.kdt == 0
    state.kD = ceil(p*state.kD)/p + 1/p;
elseif state.kD < 1 - e
    state.kD = state.kD + 1/p;
else
    state.ci = state.ci + 1;
end
end

%% ========================================================================
% Internal bracketing + fzero (so utils/ isn't required yet)
% ========================================================================
function [root, info] = dt_bracket_and_fzero(g, dt_min, dt_max, grow, C)

info = struct('dtL',NaN,'dtR',NaN,'gL',NaN,'gR',NaN,'status','');

% We want g(dtL) > 0, g(dtR) < 0 (assuming Sigma decreases with dt).
dtL = dt_min;
gL  = g(dtL);

if ~(isfinite(gL))
    root = NaN; info.status = 'g(dt_min) not finite';
    return;
end
if gL <= 0
    root = NaN; info.status = 'already <=0 at dt_min (no admissible dt)';
    info.dtL=dtL; info.gL=gL;
    return;
end

dtR = dtL;
gR  = gL;

% Expand until sign change or dt_max hit
while dtR < dt_max
    dtR = min(dt_max, dtR*grow);
    gR = g(dtR);
    if ~isfinite(gR)
        root = NaN; info.status = 'g(dt) not finite while bracketing';
        info.dtL=dtL; info.gL=gL; info.dtR=dtR; info.gR=gR;
        return;
    end
    if gR < 0
        break;
    end
    % If never crosses, loop continues
end

info.dtL = dtL; info.gL = gL;
info.dtR = dtR; info.gR = gR;

if ~(gR < 0)
    root = NaN;
    info.status = 'no sign change up to dt_max';
    return;
end

% Solve with fzero on bracket
if ~isfield(C,'opt_fzero') || isempty(C.opt_fzero)
    C.opt_fzero = optimset('Display','off');
end
root = fzero(g, [dtL, dtR], C.opt_fzero);
info.status = 'ok';
end
