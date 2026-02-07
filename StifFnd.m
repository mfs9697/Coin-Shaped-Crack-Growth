  function Stif = StifFnd(Dmat, ctx)
%STIFFND  Assemble global stiffness matrix for 3D tetrahedra (T10 here)
%         without globals.
%
% INPUT:
%   Dmat : 6x6 constitutive matrix (or 6x6 for TI plane strain style used here)
%   ctx  : struct containing mesh/assembly data:
%          ctx.ncoord, ctx.coords, ctx.nelem, ctx.eldf, ctx.nip3,
%          ctx.xip3, ctx.w3, ctx.nelnodes, ctx.connect, ctx.fixvars,
%          ctx.ndof, ctx.velem
%
% OUTPUT:
%   Stif : ndof x ndof sparse stiffness matrix with Dirichlet rows/cols pinned

ncoord   = ctx.ncoord;
coords   = ctx.coords;
nelem    = ctx.nelem;
eldf     = ctx.eldf;
nip3     = ctx.nip3;
xip3     = ctx.xip3;
w3       = ctx.w3;
nelnodes = ctx.nelnodes;
connect  = ctx.connect;
fixvars  = ctx.fixvars;
ndof     = ctx.ndof;
velem    = ctx.velem;

nfix = numel(fixvars);

m0 = ncoord * size(coords,2);
m1 = round(0.05 * m0^2);          % legacy preallocation heuristic

rw = zeros(m1,1);
cl = zeros(m1,1);
st = zeros(m1,1);
ki = 0;

% Faster lookup than "any(fixvars==rw0)" inside loops
isFixed = false(ndof,1);
isFixed(fixvars) = true;

for lmnJ = 1:nelem
    K1 = zeros(eldf);

    for jJ = 1:nip3
        B2 = B(ctx.xip3(:,jJ), lmnJ, ctx);                 % must use ctx inside B or B must access ctx another way
        K1 = K1 + w3(jJ) * (B2' * Dmat * B2) * velem(lmnJ);
    end

    for i1 = 1:nelnodes
        base_i = ncoord*(connect(i1,lmnJ)-1);

        for ii = 1:ncoord
            rw0 = base_i + ii;

            if ~isFixed(rw0)
                for j1 = 1:nelnodes
                    base_j = ncoord*(connect(j1,lmnJ)-1);

                    for jj = 1:ncoord
                        ki = ki + 1;

                        if ki > m1
                            % grow if needed (safe)
                            grow = round(0.25*m1) + 1000;
                            rw(end+grow) = 0; cl(end+grow) = 0; st(end+grow) = 0;
                            m1 = numel(rw);
                        end

                        rw(ki) = rw0;
                        cl(ki) = base_j + jj;
                        st(ki) = K1(ncoord*(i1-1)+ii, ncoord*(j1-1)+jj);
                    end
                end
            end
        end
    end
end

Stif = sparse([rw(1:ki); fixvars], ...
              [cl(1:ki); fixvars], ...
              [st(1:ki); ones(nfix,1)], ndof, ndof);
end

