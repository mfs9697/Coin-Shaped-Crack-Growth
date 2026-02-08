function cfg = cfg_Notch3D02()
%CFG_NOTCH3D02  Notch3D case configuration (parameters only).
%
% Keep it lightweight: no mesh loops, no assembly.
% This struct is intentionally compatible with Ucurr_core() expectations.

cfg = struct();

% -------------------------
% Problem meta
% -------------------------
cfg.meta.caseName = 'Notch3D02';
cfg.meta.meshFile = '';   % TODO: set your Notch mesh file path (COMSOL export) or generator output

% -------------------------
% Geometry / mesh (placeholder)
% -------------------------
% If you have a Notch mesh generator, put its knobs here.
cfg.geom = struct();  % TODO fill as needed
cfg.mesh = struct();  % TODO fill as needed

% Element / integration
cfg.elem.nelnodes = 10;     % locked to T10 in your codebase
cfg.elem.nip3     = 4;      % typical; set to what your legacy Notch used

% -------------------------
% Material (visco) parameters
% -------------------------
% ViscMod(mat) expects mat fields shown in ViscMod.m :contentReference[oaicite:5]{index=5}
cfg.mat = struct();
cfg.mat.Et_inst = 1;        % TODO
cfg.mat.Ep_inst = 1;        % TODO
cfg.mat.Gt_inst = 1;        % TODO
cfg.mat.nu32    = 0.3;      % TODO
cfg.mat.nup     = 0.3;      % TODO
cfg.mat.tau_k   = [1];      % TODO vector
cfg.mat.ratio_E3_relax = 2; % TODO
cfg.mat.ratio_E2_relax = 2; % TODO
cfg.mat.ratio_G_relax  = 2; % TODO

% -------------------------
% Cohesive law parameters (must match F0a_core usage)
% -------------------------
cfg.coh = struct();
cfg.coh.sc1 = 1;     % scaling used in F0a_core (ctx.sc1) :contentReference[oaicite:6]{index=6}
cfg.coh.a1  = 1;     % G() parameter
cfg.coh.a2  = 1;     % G() parameter
cfg.coh.Dym = 1;     % opening scale (ctx.Dym)

% -------------------------
% Control / stepping
% -------------------------
cfg.ctrl = struct();
cfg.ctrl.nt      = 100;
cfg.ctrl.p       = 5;
cfg.ctrl.dt_min  = 1e-4;
cfg.ctrl.dt_max  = 1e+1;
cfg.ctrl.dt_guess = 0.1;
cfg.ctrl.use_fzero   = true;
cfg.ctrl.bracketGrow = 2;

% Target “external” stress level for dt selection:
% In your legacy Notch snippet: sig = sc1*0.45
cfg.ctrl.sig_target = cfg.coh.sc1 * 0.45;

% -------------------------
% Solver options (inner equilibrium and outer dt)
% -------------------------
cfg.solve = struct();

% Ucurr_core currently calls:
%   fsolve(fun, [state.dU0; 0.01*state.sig], cfg.solve.fsolve_F0c) :contentReference[oaicite:7]{index=7}
cfg.solve.fsolve_F0c = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance',     1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations', 200);

cfg.solve.fsolve_dt = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance',     1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations', 50);

cfg.solve.fzero = optimset('Display','off');
end
