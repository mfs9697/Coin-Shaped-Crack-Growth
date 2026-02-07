function cfg = cfg_CoinCrack_TI_legacy()
% One place to edit inputs. No globals here.

cfg.meta.meshFile = "CoinShapedCrack_mesh_CZM.mat";

% --- Geometry / domain ---
cfg.geom.L  = 10;
cfg.geom.H  = 10;
cfg.geom.Rc = 3.0;              % = a0 (crack radius)
cfg.geom.width = 0.5 * cfg.geom.Rc;
%cfg.geom.cr_a1 = cfg.geom.Rc - 0.1*cfg.geom.width;   % geometry radii
cfg.geom.cr_a2 = cfg.geom.Rc + cfg.geom.width;

% --- Mesh controls (builder only) ---
cfg.mesh.caseType = "CZM";
cfg.mesh.nprecoh  = 0;
cfg.mesh.ncoh     = 8;
cfg.mesh.hgrad_out= 1.3;

% --- Element choice (LOCKED) ---
cfg.elem.nelnodes = 10;  % quadratic tetra only
cfg.elem.nip3     = 4;   % consistent with your current code

% --- Cohesive law (TSL) ---
cfg.czm.phi     = 1200;
cfg.czm.sig_max = 0.02;
cfg.czm.tsl_a1  = 0.002;     % TSL shape params (dimensionless in your G)
cfg.czm.tsl_a2  = 0.002;

% --- Viscoelastic / TI material ---
cfg.mat.Et_inst = 4;
cfg.mat.Ep_inst = 4;
cfg.mat.nu32    = 0.3;
cfg.mat.nup     = 0.3;
cfg.mat.Gt_inst = cfg.mat.Et_inst/2/(1+cfg.mat.nu32);

cfg.mat.tau_k = 50;
cfg.mat.ratio_E3_relax = 2;
cfg.mat.ratio_E2_relax = 2;
cfg.mat.ratio_G_relax  = 2;


% --- Solver / time loop ---
cfg.solve.fact = 100;
cfg.solve.nt   = 100;
cfg.solve.sig  = .45*cfg.czm.sig_max;
cfg.solve.p    = 5;
cfg.solve.dt0  = 1;
cfg.solve.dt_eps  = 1;
cfg.solve.Tmax    = 320;
cfg.solve.dt_max  = cfg.solve.Tmax;   % or 5*cfg.mat.tau_k(1)
cfg.solve.Nt_curve = []; % not used here, kept for later if you want curves
cfg.solve.Nscan = 3;


cfg.solve.fsolve_F0c = optimset('Display','off','Jacobian','on','MaxFunEvals',50);
cfg.solve.fsolve_Ucurr = optimoptions('fsolve','Display','iter','TolFun',1e-18,'MaxIterations',30);
% =========================================================
% Root-finding (fzero) options for dt-determination
% =========================================================
cfg.solve.fzero_dt = optimset( ...
    'Display',      'iter', ...     % 'off' | 'iter' | 'final'
    'TolX',         1e-10,  ...     % dt precision
    'MaxIter',      80,     ...     % max iterations
    'MaxFunEvals',  200 );          % max function evaluations

end
