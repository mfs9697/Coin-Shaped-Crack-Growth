function [M,Eo,Em,rhom,tEi] = ViscMod(mat)
%VISCMOD  Build TI viscoelastic moduli data from cfg.mat
%
% INPUT (mat struct expected fields):
%   mat.Et_inst           (E30)
%   mat.Ep_inst           (E10)
%   mat.Gt_inst           (G0)
%   mat.nu32
%   mat.nup               (nu21)
%   mat.tau_k             (vector), uses first value as ro
%   mat.ratio_E3_relax    (= E30 / E3i)
%   mat.ratio_E2_relax    (= E20 / E2i)
%   mat.ratio_G_relax     (= G0  / Gi)
%
% OUTPUT: same as your legacy ViscMod
%   M, Eo, Em, rhom, tEi

% --- required inputs ---
ro  = mat.ro;

% Build instantaneous moduli consistently with your earlier mapping
E30 = mat.E30;                     % instantaneous vertical
E20 = mat.E20;                     % in-plane instantaneous
G0  = mat.G0;                      % shear instantaneous

% Relaxed (branch) moduli
E3i = mat.E3i;  E31 = E30 - E3i;
E2i = mat.E2i;  E21 = E20 - E2i;
Gi  = mat.Gi;   G1  = G0  - Gi;

nu21 = mat.nu21;
nu32 = mat.nu32;

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

cc  = hm*g*(E2i*E30 - E20*E3i)^2 / aa0 / aa1 / ai;

E = [ bi/ai             bb1/aa1             -cc
      di/ai             d1/aa1              -cc
      nu32*E2i*E3i/ai   nu32*E21*E31/aa1    -cc*2*nu32
      hm*E3i^2/ai       hm*E31^2/aa1        -cc*4*g
      Gi                G1                   0
      E2i/hp/2          E21/hp/2             0 ];

ro2 = aa0/ai * ro;

M = 2;                                  % fixed in your model
rhom = zeros(6,6,M);
Em   = zeros(6,6,M);

rhom(:,:,1) = ro  * eye(6);
rhom(:,:,2) = ro2 * eye(6);

for m1 = 2:M+1
    Em(:,:,m1-1) = ...
        [E(1,m1) E(2,m1) E(3,m1) 0 0 0
         E(2,m1) E(1,m1) E(3,m1) 0 0 0
         E(3,m1) E(3,m1) E(4,m1) 0 0 0
         0 0 0 E(6,m1) 0 0
         0 0 0 0 E(5,m1) 0
         0 0 0 0 0 E(5,m1)];
end

Eonz = sum(E,2);
Eo = [Eonz(1) Eonz(2) Eonz(3) 0 0 0
      Eonz(2) Eonz(1) Eonz(3) 0 0 0
      Eonz(3) Eonz(3) Eonz(4) 0 0 0
      0 0 0 Eonz(6) 0 0
      0 0 0 0 Eonz(5) 0
      0 0 0 0 0 Eonz(5)];

tEi = [E(1,1) E(2,1) E(3,1) 0 0 0
       E(2,1) E(1,1) E(3,1) 0 0 0
       E(3,1) E(3,1) E(4,1) 0 0 0
       0 0 0 E(6,1) 0 0
       0 0 0 0 E(5,1) 0
       0 0 0 0 0 E(5,1)];
end
