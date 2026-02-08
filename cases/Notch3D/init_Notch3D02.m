function [ctx, state] = init_Notch3D02(cfg)
%INIT_NOTCH3D02  Initialize Notch3D ctx/state (globals-free path).

ctx = struct();
state = struct();

% -------------------------
% 1) Mesh / FE operators
% -------------------------
% TODO: You must load/build your Notch mesh here.
% You likely have a COMSOL mesh file, like in init_problem() :contentReference[oaicite:9]{index=9}
%
% Suggested pattern:
%   [coords, connect4, model] = Crack3D_loadMesh_COMSOL(cfg.meta.meshFile);
%   [coords, connect, ~] = T4toT10(coords, connect4);  % if needed
%   ... compute velem, quadrature, fixvars, etc.
%
% For now we error clearly if missing:
if isempty(cfg.meta.meshFile)
    error('init_Notch3D02: cfg.meta.meshFile is empty. Set it to your Notch COMSOL mesh path.');
end

% -------------------------
% 2) Viscoelastic material (ctx.M, ctx.Em, ctx.rhom, ctx.tEi)
% -------------------------
[ctx.M, ctx.Eo, ctx.Em, ctx.rhom, ctx.tEi] = ViscMod(cfg.mat);

% -------------------------
% 3) Cohesive parameters copied into ctx
% -------------------------
ctx.sc1 = cfg.coh.sc1;
ctx.a1  = cfg.coh.a1;
ctx.a2  = cfg.coh.a2;
ctx.Dym = cfg.coh.Dym;

% -------------------------
% 4) State arrays (allocate after ndof/nelem/nip3 known)
% -------------------------
% TODO after you load mesh:
% ctx.ndof, ctx.nelem, ctx.nip3, ctx.eldf must exist
% then allocate:
% state.U      = zeros(ctx.ndof,1);
% state.dU0    = zeros(ctx.ndof,1);
% state.Fpsig  = zeros(ctx.ndof,1);
% state.Ftsig  = zeros(ctx.ndof,1);
% state.S      = zeros(6,6,ctx.nip3*ctx.nelem,ctx.M);
% state.sigzip = zeros(ctx.nelem*ctx.nip3,1);
% state.t1     = zeros(cfg.ctrl.nt+1,1);

% Control indices
state.ci = 1;
state.kD = 0;

% Needed by Ucurr_core for initial guess:
state.sig = cfg.ctrl.sig_target;

end
