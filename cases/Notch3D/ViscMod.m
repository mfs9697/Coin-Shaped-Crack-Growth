function [M, Eo, Em, rhom, tEi] = ViscMod(mat)
%VISCMOD  TI viscoelastic parameters (2-term Prony-like form).
%
% Usage:
%   [M,Eo,Em,rhom,tEi] = ViscMod();        % legacy defaults
%   [M,Eo,Em,rhom,tEi] = ViscMod(mat);     % use mat fields
%
% Required/used fields in mat (defaults provided if missing):
%   ro   : base relaxation time (scalar, >0)
%   E20, E2i : transverse modulus (instant & long-time parts)
%   E30, E3i : longitudinal modulus (instant & long-time parts)
%   nu21, nu32 : Poisson-like constants
%   G_mode : 'from_E_nu' (default) or 'explicit'
%   (if explicit) G0, Gi : shear moduli (instant & long-time parts)
%
% Outputs are identical in meaning to your legacy code:
%   M    : number of relaxation branches (here fixed to 2)
%   Eo   : instantaneous stiffness matrix (6x6)
%   Em   : tangent moduli blocks per branch (6x6xM)
%   rhom : relaxation times arranged as 6x6xM (legacy pattern)
%   tEi  : long-time stiffness matrix (6x6)

if nargin < 1 || isempty(mat)
    mat = struct();
end

% -------------------------
% Defaults (match your current hard-coded set)
% -------------------------
mat = set_default(mat, 'ro',   20);

mat = set_default(mat, 'E20',  4);
mat = set_default(mat, 'E2i',  1);

mat = set_default(mat, 'E30',  4);
mat = set_default(mat, 'E3i',  1);

mat = set_default(mat, 'nu21', 0.3);
mat = set_default(mat, 'nu32', 0.3);

mat = set_default(mat, 'G_mode', 'from_E_nu');  % or 'explicit'

% -------------------------
% Unpack (legacy names)
% -------------------------
ro   = mat.ro;

E20  = mat.E20;  E2i = mat.E2i;  E21 = E20 - E2i;
E30  = mat.E30;  E3i = mat.E3i;  E31 = E30 - E3i;

nu21 = mat.nu21;
nu32 = mat.nu32;

G0 = mat.G0;
Gi = mat.Gi;
G1 = G0 - Gi;

% -------------------------
% Sanity checks
% -------------------------
must_be_pos_scalar(ro, 'ro');

% Basic physical plausibility checks (not too strict, but catches typos)
if any([E20,E2i,E30,E3i] <= 0)
    error('ViscMod:BadE', 'E20,E2i,E30,E3i must be > 0.');
end
if E2i > E20 || E3i > E30
    error('ViscMod:BadSplit', 'Long-time parts E2i<=E20 and E3i<=E30 required.');
end

% -------------------------
% Legacy algebra (UNCHANGED)
% -------------------------
g  = nu32^2;
hp = 1 + nu21;
hm = 1 - nu21;

aa0 = E30*hm - 2*E20*g;
aa1 = E31*hm - 2*E21*g;

bb1 = E21*(E31 - E21*g)/hp;
d1  = E21*(E31*nu21 + E21*g)/hp;

ai  = E3i*hm - 2*E2i*g;
bi  = E2i*(E3i - E2i*g)/hp;
di  = E2i*(E3i*nu21 + E2i*g)/hp;

cc = hm*g*(E2i*E30 - E20*E3i)^2/aa0/aa1/ai;

E = [ bi/ai               bb1/aa1             -cc
      di/ai               d1/aa1              -cc
      nu32*E2i*E3i/ai     nu32*E21*E31/aa1    -cc*2*nu32
      hm*E3i^2/ai         hm*E31^2/aa1        -cc*4*g
      Gi                  G1                  0
      E2i/hp/2            E21/hp/2            0 ];

ro2 = aa0/ai * ro;

M = 2;

% -------------------------
% rhom and Em (legacy pattern)
% -------------------------
rhom = zeros(6,6,M);
Em   = zeros(6,6,M);

rhom(:,:,1) = [ro ro ro 0  0  0
               ro ro ro 0  0  0
               ro ro ro 0  0  0
               0  0  0  ro 0  0
               0  0  0  0  ro 0
               0  0  0  0  0  ro];

rhom(:,:,2) = [ro2 ro2 ro2 0 0 0
               ro2 ro2 ro2 0 0 0
               ro2 ro2 ro2 0 0 0
               0   0   0   0 0 0
               0   0   0   0 0 0
               0   0   0   0 0 0];

for m = 2:(M+1)
    Em(:,:,m-1) = [E(1,m) E(2,m) E(3,m) 0      0      0
                   E(2,m) E(1,m) E(3,m) 0      0      0
                   E(3,m) E(3,m) E(4,m) 0      0      0
                   0      0      0      E(6,m) 0      0
                   0      0      0      0      E(5,m) 0
                   0      0      0      0      0      E(5,m)];
end

% Instantaneous moduli
Eonz = sum(E,2);
Eo = [Eonz(1) Eonz(2) Eonz(3) 0       0       0
      Eonz(2) Eonz(1) Eonz(3) 0       0       0
      Eonz(3) Eonz(3) Eonz(4) 0       0       0
      0       0       0       Eonz(6) 0       0
      0       0       0       0       Eonz(5) 0
      0       0       0       0       0       Eonz(5)];

% Long-time moduli
tEi = [E(1,1) E(2,1) E(3,1) 0      0      0
       E(2,1) E(1,1) E(3,1) 0      0      0
       E(3,1) E(3,1) E(4,1) 0      0      0
       0      0      0      E(6,1) 0      0
       0      0      0      0      E(5,1) 0
       0      0      0      0      0      E(5,1)];

% Final sanity
if ~all(isfinite(Eo(:))) || ~all(isfinite(tEi(:))) || ~all(isfinite(Em(:))) || ~all(isfinite(rhom(:)))
    error('ViscMod:NonFinite', 'Non-finite values produced; check inputs.');
end
if any(rhom(:) < 0)
    error('ViscMod:BadRho', 'Negative relaxation times produced (rhom).');
end

end

% ===== helpers =====

function s = set_default(s, name, val)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = val;
end
end

function must_be_pos_scalar(x, name)
if ~isscalar(x) || ~isfinite(x) || x <= 0
    error('ViscMod:BadScalar', '%s must be a positive finite scalar.', name);
end
end
