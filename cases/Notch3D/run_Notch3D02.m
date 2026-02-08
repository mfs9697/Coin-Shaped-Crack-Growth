function run_Notch3D02()
%RUN_NOTCH3D02  Step2-style driver skeleton for Notch3D case.
%
% Goals:
%   1) Make the workflow explicit: cfg -> init -> time loop -> dt selection -> update -> plots
%   2) Keep physics/numerics unchanged (initially)
%   3) Prepare the codebase for replacing dt-fsolve by bracketed fzero later
%
% Expected companion files (recommended):
%   cases/Notch3D/cfg_Notch3D02.m
%   cases/Notch3D/init_Notch3D02.m
% Optional (later):
%   utils/dt_bracket_fzero.m

clc;

%% ------------------------------------------------------------
% 0) Paths (keep minimal; replace by startup.m later if you want)
% -------------------------------------------------------------
repoRoot = fileparts(fileparts(mfilename('fullpath'))); % .../cases
projRoot = fileparts(repoRoot);                         % project root
addpath(genpath(projRoot));

%% -------------------------
% 1) Configuration (cfg)
% --------------------------
C = cfg_Notch3D02();
[ctx, state, C] = init_Notch3D02(C);

% Make output directory (optional but recommended)
if ~isfield(C,'outDir') || isempty(C.outDir)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    C.outDir = fullfile(projRoot, 'out', ['Notch3D02_' ts]);
end
if ~exist(C.outDir,'dir'); mkdir(C.outDir); end

logFile = fullfile(C.outDir, 'log.txt');
fid = fopen(logFile,'w');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'run_Notch3D02 started: %s\n', datestr(now));
fprintf(fid, 'MATLAB version: %s\n\n', version);

%% -------------------------
% 2) Initialization (init)
% --------------------------
% init should build:
%   ctx: constants (mesh, operators, material constants, indices)
%   state: evolving variables (U, S, Fpsig, Ftsig, sigzip, ci, kD, ...)
[ctx, state] = init_Notch3D02(C);

% Optional sanity prints
fprintf(fid, 'ndof=%d, nelem=%d, nnodes=%d\n', ctx.ndof, ctx.nelem, ctx.nnodes);
fprintf(fid, 'Initial ci=%d, kD=%.6g\n\n', state.ci, state.kD);

%% -----------------------------------------
% 3) Main loop (incremental quasi-static)
% ------------------------------------------
for kdt = 0:C.nt

    state.kdt = kdt;

    % -----------------------------------------
    % 3a) Initial step (kdt == 0)
    % -----------------------------------------
    if kdt == 0
        % === Put your existing "elastic initialization" block here ===
        %
        % Example mapping from your legacy code:
        %   K = Stif(Eo);
        %   dU = fsolve(@F0c, dU0, opt);  dU = full(dU);
        %   kD = 2*dU(vc(1))/Dym;
        %
        % In this refactor:
        %   - store any "K" into ctx (or state) as needed
        %   - store dU into state.dU
        %   - store kD into state.kD
        %   - keep dU0 updated in state.dU0 for future inner solves

        [ctx, state] = do_initial_elastic_step(ctx, state, C, fid);

    else
        % -----------------------------------------
        % 3b) dt selection for kdt > 0
        % -----------------------------------------
        % Legacy behavior:
        %   dt0  = 1e-4;
        %   sig0 = Ucurr(dt0,kD,ci);
        %   if sig0 > sig, dt = fsolve(@(dt)Ucurr(dt,kD,ci)-sig, 0.1, opt2); else break; end
        %
        % Step2-style: isolate the scalar mismatch g(dt) = Sigma(dt)-sig_target

        sig_target = C.sig_target;   % could be constant or schedule
        dt_min     = C.dt_min;
        dt_max     = C.dt_max;

        % --- quick probe at dt_min ---
        sig_min = Sigma_of_dt(dt_min, ctx, state, C);   % calls your Ucurr-equivalent
        fprintf(fid, 'kdt=%d: Sigma(dt_min)=%.6e, sig_target=%.6e\n', kdt, sig_min, sig_target);

        if sig_min <= sig_target
            fprintf(fid, 'Stop: Sigma(dt_min) <= target (no admissible quasi-static dt under this control)\n');
            break;
        end

        % --- dt solve (Stage 1: keep legacy fsolve; Stage 2: replace by bracketing+fzero) ---
        if isfield(C,'use_fzero') && C.use_fzero
            % Preferred (later): bracket + fzero
            g = @(dt) Sigma_of_dt(dt, ctx, state, C) - sig_target;
            [dt, info] = dt_bracket_fzero(g, dt_min, dt_max, C.bracketGrow, C.opt_fzero);

            fprintf(fid,'  bracket: dtL=%.3e gL=%.3e | dtR=%.3e gR=%.3e | status=%s\n',...
                info.dtL, info.gL, info.dtR, info.gR, info.status);

            if ~isfinite(dt)
                fprintf(fid,'Stop: no root bracketed in [dt_min, dt_max].\n');
                break;
            end
        else
            % Legacy: scalar fsolve on dt (kept temporarily to preserve behavior)
            opt2 = C.opt_dt_fsolve;
            dt_guess = C.dt_guess;
            dt = fsolve(@(dt) Sigma_of_dt(dt, ctx, state, C) - sig_target, dt_guess, opt2);
        end

        state.dt = dt;

        % -----------------------------------------
        % 3c) Advance state (compute dU, update U, update visco history)
        % -----------------------------------------
        % This should incorporate your existing big block:
        %   U = U + dU
        %   update S at integration points
        %   compute nodal stresses (sigz)
        %   update Fpsig/Ftsig
        %   update time vector t1
        %
        % Keep exactly the same numerical operations initially.

        [ctx, state, diag] = advance_one_increment(ctx, state, C, fid);

        % -----------------------------------------
        % 3d) Plot / output (optional but useful)
        % -----------------------------------------
        if C.doPlot
            do_plots(ctx, state, C);
        end

        % -----------------------------------------
        % 3e) Update (ci, kD) stepping rule
        % -----------------------------------------
        % Keep your exact rule:
        %   if kdt==0, kD=ceil(p*kD)/p+1/p;
        %   elseif kD<1-e, kD=kD+1/p;
        %   else, ci=ci+1;
        %   end
        state = update_control_indices(ctx, state, C);

    end

    % optional pause
    if isfield(C,'pauseSec') && C.pauseSec > 0
        pause(C.pauseSec);
    end
end

fprintf(fid, '\nrun_Notch3D02 finished: %s\n', datestr(now));
end

%% ========================================================================
% Helper stubs (to be filled with code moved from Notch3D02_legacy.m)
% ========================================================================

function [ctx, state] = do_initial_elastic_step(ctx, state, C, fid)
% Move your kdt==0 block here, unchanged.
% Must set:
%   state.dU, state.kD, ctx.K (or wherever K is stored), state.dU0
error('do_initial_elastic_step: TODO (move legacy kdt==0 block here).');
end

function sig = Sigma_of_dt(dt, ctx, state, cfg)
    [sig, state2] = Ucurr_core(dt, cfg, ctx, state);
    % IMPORTANT: Sigma_of_dt should NOT mutate the caller state during probing,
    % otherwise dt bracketing becomes path-dependent. So we discard state2 here.
end


function [ctx, state, diag] = advance_one_increment(ctx, state, C, fid)
% Move the big per-element loop and all updates here, unchanged.
% Must update:
%   state.U, state.sigzip, state.S, state.Fpsig, state.Ftsig, state.t1, ...
diag = struct();
error('advance_one_increment: TODO (move legacy increment update block here).');
end

function do_plots(ctx, state, C)
% Move your figure(1)/figure(3) code here (optional).
end

function state = update_control_indices(ctx, state, C)
% Implement your kD/ci stepping rule here, unchanged initially.
% Needs C.p and ctx.e (tolerance) possibly.
error('update_control_indices: TODO (move legacy kD/ci update rule here).');
end
