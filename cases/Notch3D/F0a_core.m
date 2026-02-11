function [f, Jac] = F0a_core(us, ctx, state)
%F0A_CORE  Coupled residual/Jacobian for [du; sig1] without globals.
%
% us    : (ndof+1) x 1  -> [du; sig1]
% ctx   : constant data for this solve (mesh/assembly/CZM params)
% state : current state (U, Fpsig, Ftsig, crack-front index, kD, etc.)

ndof = ctx.ndof;
du   = us(1:ndof);
sig1 = us(end);

% --- unpack (explicit) ---
K        = ctx.K;

nfnodes  = ctx.nfnodes;
nfnodes2 = ctx.nfnodes2;
me       = ctx.me;
mc       = ctx.mc;

je       = ctx.je;
jc       = ctx.jc;
re       = ctx.re;

cfar     = ctx.cfar;          % cohesive face area
Rface    = ctx.R;             % face operator
indq0    = ctx.indq0;

Dym      = ctx.Dym;
sc1      = ctx.sc1;
a1       = ctx.a1;
a2       = ctx.a2;

vc = state.vc;

U        = state.U;
ci       = state.ci;
kD       = state.kD;
Fpsig    = state.Fpsig;
Ftsig    = state.Ftsig;

% --- cohesive opening / traction ---
D1   = 2*(U(jc) + du(jc)) / Dym;      % dimensionless opening measure
tmp1 = -sc1 * G(D1, a1, a2);

t  = tmp1(:,1);
ta = tmp1(:,2) * 2 / Dym;             % chain rule

% --- cohesive traction vector and cohesive tangent block ---
rc  = zeros(nfnodes*mc, 1);
p1  = zeros(nfnodes2*mc, 1);
q1  = p1;
jv  = p1;                              % values for sparse() (rename to avoid shadowing Jac)

for i4 = 1:mc
    ind3 = nfnodes*(i4-1) + (1:nfnodes);

    R1 = cfar(i4) * Rface;
    rc(ind3) = R1 * t(ind3);

    indp = nfnodes2*(i4-1) + (1:nfnodes2);
    indq = nfnodes2*(i4-1) + indq0;

    p1(indp) = repmat(jc(ind3), nfnodes, 1);
    q1(indq) = p1(indp);

    jv(indp) = -reshape(R1 .* repmat(ta(ind3)', nfnodes, 1), nfnodes2, 1);
end

% --- global RHS vector (external + cohesive) ---
r = sparse([je; jc], ones(nfnodes*(me+mc),1), [sig1*re; rc], ndof, 1);

% --- residual (ndof + 1 equations) ---
f = [ K*du - (r - Fpsig + Ftsig) ; ...
      U(vc(ci)) + du(vc(ci)) - kD*Dym/2 ];

% --- Jacobian (ndof+1) x (ndof+1) ---
[rw, cl, st] = find(K);

Jac = sparse( ...
    [rw;   p1;  je;               ndof+1], ...
    [cl;   q1;  (ndof+1)*ones(me*nfnodes,1); vc(ci)], ...
    [st;   jv;  -re;              1], ...
    ndof+1, ndof+1 );
end
