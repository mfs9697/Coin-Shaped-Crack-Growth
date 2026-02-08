function C = cfg_Notch3D02()
%CFG_NOTCH3D02  Configuration for Notch3D case (parameters only).
%
% Keep this file free of heavy computations and mesh loops.
% It should only define numeric parameters, switches, and solver options.

C = struct();

% --- time stepping / continuation ---
C.nt          = 100;     % number of increments
C.dt_min      = 1e-4;    % strictly positive
C.dt_max      = 1e+1;    % cap for bracketing (tune later)
C.dt_guess    = 0.1;     % legacy (used only if you keep dt-fsolve temporarily)
C.bracketGrow = 2;       % dt growth factor for bracketing
C.use_fzero   = true;    % switch outer dt solver: true=fzero, false=fsolve

% --- control rule ---
C.p      = 5;            % kD stepping granularity (legacy)
C.epsTol = 1e-12;        % small tolerance used in comparisons

% --- load level (target stress / load factor) ---
% In your legacy code: sig = sc1 * .45;  (if sc1 is known only after init, set later)
C.sig_factor = 0.45;     % multiply by sc1 to get target sigma

% --- plotting / output ---
C.doPlot   = true;
C.pauseSec = 0.0;

% --- solver options (keep close to legacy) ---
C.opt_inner = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',200);

C.opt_dt_fsolve = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations',50);

C.opt_fzero = optimset('Display','off');  % fzero uses optimset
end
