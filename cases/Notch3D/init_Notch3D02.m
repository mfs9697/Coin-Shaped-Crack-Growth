function [ctx, state, C] = init_Notch3D02(C)
%INIT_NOTCH3D02  Build ctx/state for Notch3D case.
%
% ctx   : constants (mesh, operators, material parameters, indices)
% state : evolving vars (U, S, Fpsig, Ftsig, ci, kD, etc.)
%
% IMPORTANT: Initially, it’s okay if you keep some legacy globals,
% but the goal is to progressively move them into ctx/state.

ctx   = struct();
state = struct();

% -------------------------
% 1) Load/build mesh & FE data
% -------------------------
% TODO: move your legacy mesh import/build here:
% - coords, connect, nnodes, nelem, ndof, ...
% - nip3, xip3, w3, velem, BN handle, Nextr, etc.
%
% Example placeholders:
% ctx.coords   = coords;
% ctx.connect  = connect;
% ctx.nnodes   = nnodes;
% ctx.nelem    = nelem;
% ctx.ndof     = ndof;

% -------------------------
% 2) Cohesive index mapping (ci, vc, jc, etc.)
% -------------------------
% TODO: move your legacy cohesive mapping here:
% - compute cohlsym/cohlsym sorting, vcls, etc.
% - set vc, jc, etc.

% -------------------------
% 3) Viscoelastic material precompute
% -------------------------
% Legacy block:
% [M,Eo,Em,rhom,tEi]=ViscMod();
% zetamn=rhom; etamn=rhom;
%
% Store them in ctx or state:
% ctx.M = M; ctx.Eo = Eo; ctx.Em = Em; ctx.rhom = rhom; ctx.tEi = tEi;

% -------------------------
% 4) Initialize evolving state arrays
% -------------------------
% Legacy:
% U=zeros(ndof,1); dU=U; Fpsig=U; Ftsig=U;
% S=zeros(6,6,nip3*nelem,M);
% sigzip=zeros(nelem*nip3,1); ...
%
% state.U      = zeros(ctx.ndof,1);
% state.dU0    = zeros(ctx.ndof,1);      % initial guess for inner solve
% state.Fpsig  = zeros(ctx.ndof,1);
% state.Ftsig  = zeros(ctx.ndof,1);
% state.S      = zeros(6,6,ctx.nip3*ctx.nelem,ctx.M);
% state.sigzip = zeros(ctx.nelem*ctx.nip3,1);
% state.t1     = zeros(C.nt+1,1);

% -------------------------
% 5) Initial control values
% -------------------------
state.ci = 1;
state.kD = NaN;  % will be set after initial elastic step

% -------------------------
% 6) Define target load now if it depends on sc1
% -------------------------
% If sc1 is known only after mesh/material init:
% C.sig_target = ctx.sc1 * C.sig_factor;
% else set directly in cfg.

end
