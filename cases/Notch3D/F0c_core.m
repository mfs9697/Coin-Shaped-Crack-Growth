function [f, Jac] = F0c_core(du, ctx, state)
%F0C_CORE  First-step residual and Jacobian without globals.
%
% du    : ndof x 1 increment
% ctx   : struct with geometry/assembly data (constant in time step)
% state : struct with current state (U, external load level, etc.)
%
% Returns:
%   f   : ndof x 1
%   Jac : ndof x ndof (sparse)

% --- unpack what we need (explicitly) ---
ndof     = ctx.ndof;
K        = ctx.K;
U        = state.U;

% cohesive mapping / quadrature for faces
nfnodes  = ctx.nfnodes;
nfnodes2 = ctx.nfnodes2;
mc       = ctx.mc;

jc       = ctx.jc;      % dof indices on cohesive faces (stacked)
cfar     = ctx.cfar;    % cohesive face areas (your note)
Rface    = ctx.R;       % face “mass” matrix
indq0    = ctx.indq0;

% external load mapping
je       = ctx.je;
re       = ctx.re;
me       = ctx.me;
sig_ext  = state.sig;   % applied traction level (scalar)

% cohesive law scaling
Dym      = ctx.Dym;
sc1      = ctx.sc1;
a1       = ctx.a1;
a2       = ctx.a2;

% optional history forces (kept explicit)
Fpsig    = state.Fpsig;
Ftsig    = state.Ftsig;

% --- cohesive opening measure on cohesive DOFs ---
% (same expression you use in F0a)
D1   = 2*(U(jc) + du(jc)) / Dym;         % dimensionless opening measure
tmp1 = -sc1 * G(D1, a1, a2);             % bT(:,1)=g, bT(:,2)=dg/dy

t  = tmp1(:,1);                          % traction shape (dimensionless)
ta = tmp1(:,2) * 2 / Dym;                % dg/d(physical opening) chain rule

% --- assemble cohesive nodal tractions rc and cohesive tangent ---
rc  = zeros(nfnodes*mc,1);
p1  = zeros(nfnodes2*mc,1);
q1  = p1;
jv  = p1;                                 % values for sparse Jacobian block

for i4 = 1:mc
    ind3 = nfnodes*(i4-1) + (1:nfnodes);

    % face operator scaled by cohesive face area
    R1 = cfar(i4) * Rface;

    % nodal cohesive traction vector contribution
    rc(ind3) = R1 * t(ind3);

    % Jacobian assembly indices (same pattern as your F0a)
    indp = nfnodes2*(i4-1) + (1:nfnodes2);
    indq = nfnodes2*(i4-1) + indq0;

    % row/col indices in global dof space
    p1(indp) = repmat(jc(ind3), nfnodes, 1);
    q1(indq) = p1(indp);

    % values: d(rc)/d(du) = -R1 .* ta  (consistent with your sign convention)
    jv(indp) = -reshape(R1 .* repmat(ta(ind3)', nfnodes, 1), nfnodes2, 1);
end

% --- assemble external + cohesive RHS into global ndof vector ---
r = sparse([je; jc], ones(nfnodes*(me+mc),1), [sig_ext*re; rc], ndof, 1);

% --- residual ---
% same structure you use later: K*du - (r - Fpsig + Ftsig)
f = K*du - (r - Fpsig + Ftsig);

% --- Jacobian ---
% start from K, add cohesive tangent block
[rw, cl, st] = find(K);
Jac = sparse([rw; p1], [cl; q1], [st; jv], ndof, ndof);
end
