function Bout = B(xi, lmn, ctx)
%B  Strain-displacement matrix for 3D tetrahedral element (T10 supported).
%
% Bout = B(xi, lmn, ctx)
%   xi  : 4x1 barycentric coordinates [x1;x2;x3;x4]
%   lmn : element index
%   ctx : struct with fields:
%         ctx.nelnodes, ctx.ncoord, ctx.coords, ctx.connect, ctx.eldf

nelnodes = ctx.nelnodes;
ncoord   = ctx.ncoord;
coords   = ctx.coords;
connect  = ctx.connect;
eldf     = ctx.eldf;

x1 = xi(1); x2 = xi(2); x3 = xi(3); x4 = xi(4);

dNdxi = zeros(nelnodes, ncoord);

if nelnodes == 10
    dNdxi(1,1)        = x1 - 0.25;
    dNdxi(2,2)        = x2 - 0.25;
    dNdxi(3,3)        = x3 - 0.25;
    dNdxi(4,1:3)      = -x4 + 0.25;

    dNdxi(5,[1,2])    = [x2, x1];
    dNdxi(6,[2,3])    = [x3, x2];
    dNdxi(7,[3,1])    = [x1, x3];

    dNdxi(8,[1,2,3])  = [x4-x1, -x1,   -x1];
    dNdxi(9,[2,3,1])  = [x4-x2, -x2,   -x2];
    dNdxi(10,[3,1,2]) = [x4-x3, -x3,   -x3];

    dNdxi = 4*dNdxi;
else
    error('B(): only nelnodes==10 (quadratic tetra) is supported in this refactor.');
end

% Jacobian mapping: x = sum N_i * x_i  => dx/dxi = X * dN/dxi
X = coords(:, connect(:,lmn));          % 3 x nelnodes
dxdxi = X * dNdxi;                      % 3 x 3

% Safer than inv()
dNdx = dNdxi / dxdxi;                   % nelnodes x 3

Bout = zeros(6, eldf);

for i1 = 1:nelnodes
    b = dNdx(i1,:);  % [dN/dx dN/dy dN/dz]
    col = ncoord*(i1-1) + (1:ncoord);

    Bout(:, col) = [ b(1) 0    0
                     0    b(2) 0
                     0    0    b(3)
                     b(2) b(1) 0
                     b(3) 0    b(1)
                     0    b(3) b(2) ];
end
end
