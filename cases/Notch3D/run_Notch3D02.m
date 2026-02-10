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
[ctx, state, C] = init_Notch3D02(C);

%% Output folder + log
if ~isfield(C,'outDir') || isempty(C.outDir)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    C.outDir = fullfile(projRoot, 'out', ['Notch3D02_' ts]);
end
if ~exist(C.outDir,'dir'); mkdir(C.outDir); end

logFile = fullfile(C.outDir, 'log.txt');
fid = fopen(logFile,'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'run_Notch3D02 started: %s\n', datestr(now));
fprintf(fid, 'MATLAB version: %s\n\n', version);

fprintf(fid, 'ndof=%d, nelem=%d, nnodes=%d, nip3=%d, M=%d\n', ...
    ctx.ndof, ctx.nelem, ctx.nnodes, ctx.nip3, ctx.M);
fprintf(fid, 'Initial ci=%d\n\n', state.ci);

%% Main loop
for kdt = 0:C.nt
    state.kdt = kdt;

    if kdt == 0
        % Elastic init (matches legacy: K(Eo), then fsolve(F0c), then kD)
        [ctx, state] = do_initial_elastic_step(ctx, state, C, fid);
        % After kdt==0, update stepping rule once (legacy did this after plotting)
        state = update_control_indices(ctx, state, C);
        continue;
    end

    % ---- dt selection (Sigma(dt)=sig_target) ----
    sig_target = C.sig_target;

    dt_min = C.dt_min;
    dt_max = C.dt_max;

    sig_min = Sigma_of_dt(dt_min, ctx, state, C);
    fprintf(fid, 'kdt=%d: Sigma(dt_min)=%.6e, target=%.6e\n', kdt, sig_min, sig_target);

    % If Sigma(dt_min) <= target, no admissible quasi-static dt (under this control)
    if sig_min <= sig_target
        fprintf(fid, 'Stop: Sigma(dt_min) <= target (no admissible quasi-static dt)\n');
        break;
    end

    if isfield(C,'use_fzero') && C.use_fzero
        g = @(dt) Sigma_of_dt(dt, ctx, state, C) - sig_target;
        [dt, info] = dt_bracket_and_fzero(g, dt_min, dt_max, C.bracketGrow, C);
        fprintf(fid,'  bracket: dtL=%.3e gL=%.3e | dtR=%.3e gR=%.3e | %s\n', ...
            info.dtL, info.gL, info.dtR, info.gR, info.status);
        if ~isfinite(dt)
            fprintf(fid,'Stop: no root bracketed in [dt_min, dt_max].\n');
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
    [ctx, state, diag] = advance_one_increment(ctx, state, C, fid); %#ok<NASGU>

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

fprintf(fid, '\nrun_Notch3D02 finished: %s\n', datestr(now));
end

%% ========================================================================
% kdt==0 elastic step (legacy-faithful via F0c_core)
% ========================================================================
function [ctx, state] = do_initial_elastic_step(ctx, state, C, fid)

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

fprintf(fid,'kdt=0: sig_ext=%.6e | kD=%.6e\n', state.sig, state.kD);
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
function [ctx, state, diag] = advance_one_increment(ctx, state, C, fid)
diag = struct();

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
state.U  = state.U + dU;
state.dU0 = dU;  % next initial guess

% 3) Update internal history S at integration points (legacy scheme)
M     = ctx.M;
indU  = zeros(ctx.eldf,1);

for el = 1:ctx.nelem
    nod = ctx.connect(:,el);
    for j = 1:ctx.nelnodes
        indU((3*j-2):(3*j)) = (3*nod(j)-2):(3*nod(j));
    end

    for ip = 1:ctx.nip3
        B3 = B(ctx.xip3(:,ip), el, ctx);
        deps = B3 * dU(indU);              % 6x1
        lmnip = (el-1)*ctx.nip3 + ip;

        deps1 = repmat(deps', 6, 1);
        for m = 1:M
            state.S(:,:,lmnip,m) = state.zetamn(:,:,m).*state.S(:,:,lmnip,m) + ...
                state.etamn(:,:,m).*(1 - state.zetamn(:,:,m)).*deps1;
        end
    end
end

% 4) Update force history (legacy: Fpsig = K*dU + Fpsig - Ftsig)
state.Fpsig = state.K*dU + state.Fpsig - state.Ftsig;

% 5) Update time vector
kdt = state.kdt;
if kdt == 1
    state.t1(2) = state.dt;
else
    state.t1(kdt+1) = state.t1(kdt) + state.dt;
end

fprintf(fid,'kdt=%d: dt=%.4e | sig=%.6e | ci=%d | kD=%.6e\n', ...
    kdt, state.dt, state.sig, state.ci, state.kD);

diag.sig = state.sig;
diag.dt  = state.dt;
end

%% ========================================================================
% Plot hook (optional)
% ========================================================================
function do_plots(ctx, state, C) %#ok<INUSD>
% Keep empty for now, or move your legacy figure(1)/figure(3) plotting here.
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
