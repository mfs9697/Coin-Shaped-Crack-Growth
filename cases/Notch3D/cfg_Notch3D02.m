function cfg = cfg_Notch3D02()
%CFG_NOTCH3D02  Legacy-conform configuration for Notch3D02 case.
% One place to edit inputs. No globals here.
%
% Conforms to cfg_CoinCrack_TI_legacy style (cfg.solve.*),
% but also provides Step2-style aliases used in run_Notch3D02/init_Notch3D02.

cfg = struct();

% =========================================================
% Meta / files
% =========================================================
cfg.meta.caseName = 'Notch3D02';


% Geometry files are hardcoded in init_Notch3D02.m as:
%   cases/Notch3D/Geometry/crd05.txt
%   cases/Notch3D/Geometry/con05.txt
% so no need to specify meshFile here unless you want.

% =========================================================
% Geometry constants used in legacy BC logic
% =========================================================
cfg.geom.e    = 1e-8;
cfg.geom.b1   = 1.5;
cfg.geom.lcoh = 3.0;


% =========================================================
% Viscoelastic parameters
% =========================================================

cfg.mat.ro  = 20;
cfg.mat.E20 = 4;  cfg.mat.E2i = 1;
cfg.mat.E30 = 4;  cfg.mat.E3i = 1;
cfg.mat.nu21 = 0.3;
cfg.mat.nu32 = 0.3;

cfg.mat.G0=cfg.mat.E30/2/(1+cfg.mat.nu32); 
cfg.mat.Gi=cfg.mat.E3i/2/(1+cfg.mat.nu32); 
cfg.mat.G1=cfg.mat.G0-cfg.mat.Gi;

% =========================================================
% Cohesive law (TSL) parameters 
% =========================================================

cfg.coh = struct();

cfg.coh.phi1 = 150;

cfg.coh.a1 = 0.001;
cfg.coh.a2 = cfg.coh.a1;

cfg.coh.sc1 = 4/1000;

c = (3 - 2*cfg.coh.a1 + 3*cfg.coh.a2)/6;
cfg.coh.Dym = cfg.coh.phi1 / c / cfg.coh.sc1 * 1e-7;

% ---- flat aliases expected by init_Notch3D02 / F0a_core ----
cfg.sc1 = cfg.coh.sc1;
cfg.a1  = cfg.coh.a1;
cfg.a2  = cfg.coh.a2;
cfg.Dym = cfg.coh.Dym;


% =========================================================
% Solver / time loop (legacy-style block)
% =========================================================
cfg.solve = struct();

cfg.solve.nt  = 100;          % legacy: nt=100
cfg.solve.p   = 5;            % legacy: p=5

% legacy Notch: sig = sc1*0.45
cfg.solve.sig = 0.45 * cfg.coh.sc1;

% dt controls
cfg.solve.dt0     = 1e-4;     % legacy Notch probe dt0=1e-4
cfg.solve.dt_max  = 1e+1;     % safe cap (can be larger if needed)
cfg.solve.dt_guess = 0.1;     % legacy-ish initial guess for dt solve

% A small tolerance used in control stepping (kD < 1-e)
cfg.solve.epsTol = 1e-12;

% =========================================================
% Root-finding controls (dt via fzero)
% =========================================================
cfg.solve.use_fzero   = true;
cfg.solve.bracketGrow = 2.0;

cfg.solve.fzero_dt = optimset( ...
    'Display',      'off', ...
    'TolX',         1e-10, ...
    'MaxIter',      80, ...
    'MaxFunEvals',  200 );

% =========================================================
% Inner equilibrium solver options (F0c/F0a in your code)
% =========================================================
% In CoinCrack legacy:
%   optimset('Display','off','Jacobian','on','MaxFunEvals',50)
% Keep same spirit.
cfg.solve.fsolve_F0c = optimset( ...
    'Display',   'iter', ...
    'Jacobian',  'on', ...
    'MaxFunEvals', 80, ...
    'TolFun',    1e-12, ...
    'TolX',      1e-12 );

% Optional: if you still keep dt-fsolve temporarily
cfg.solve.fsolve_dt = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance',     1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxIterations', 50);

% =========================================================
% Step2-style aliases (so run_Notch3D02 can use one namespace)
% =========================================================
cfg.nt          = cfg.solve.nt;
cfg.p           = cfg.solve.p;

cfg.sig_target  = cfg.solve.sig;       % used in init + driver
cfg.dt_min      = cfg.solve.dt0;
cfg.dt_max      = cfg.solve.dt_max;
cfg.dt_guess    = cfg.solve.dt_guess;

cfg.use_fzero   = cfg.solve.use_fzero;
cfg.bracketGrow = cfg.solve.bracketGrow;
cfg.opt_fzero   = cfg.solve.fzero_dt;

cfg.epsTol      = cfg.solve.epsTol;

% If your driver uses these names:
cfg.opt_dt_fsolve = cfg.solve.fsolve_dt;

cfg.doPlot = 1; 
end
