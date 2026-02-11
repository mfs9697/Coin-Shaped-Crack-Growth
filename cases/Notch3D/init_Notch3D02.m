function [ctx, state] = init_Notch3D02(C)
%INIT_NOTCH3D02  Build ctx/state for Notch3D case (globals-free).
%
% Matches Notch3D02_legacy boundary conditions:
%   fixed1: z DOF on z≈0 and y>b1+lcoh
%   fixed2: y DOF on y≈0
%   fixed3: x DOF on x≈0
% Cohesive plane faces are restricted to y<b1+lcoh region.

ctx = struct();
state = struct();

% -------------------------
% 0) Paths
% -------------------------
here = fileparts(mfilename('fullpath'));
geomDir = fullfile(here,'Geometry');

crdFile = fullfile(geomDir,'crd05.txt');
conFile = fullfile(geomDir,'con05.txt');

% -------------------------
% 1) Load geometry (T4) and convert to T10
% -------------------------
coords_in   = readmatrix(crdFile);     % [nnodes x 3]
connect4_in = readmatrix(conFile);     % [nelem  x 4]

model = createpde();
geometryFromMesh(model,coords_in,connect4_in);

figure(2); clf(2); 
subplot(1,2,1); 
pdegplot(model,"FaceAlpha",0.9); hold on
view([120,30])
subplot(1,2,2); 
pdemesh(model,"FaceAlpha",0.9)
view([120,30])

[coords, connect] = T4toT10(coords_in, connect4_in);   % 10 x nelem

ctx.coords   = coords;
ctx.connect  = connect;
ctx.nnodes   = size(coords,2);
ctx.nelem    = size(connect,2);
ctx.nelnodes = 10;
ctx.ndof     = 3*ctx.nnodes;
ctx.eldf     = 3*ctx.nelnodes;
ctx.nip3     = 4;

% -------------------------
% 2) Integration points / weights
% -------------------------
[ctx.xip3, ctx.w3, ctx.Nextr] = IntegrationParameters(ctx.nip3, ctx.nelnodes, @N_3DT);

% -------------------------
% 3) Element volumes
% -------------------------
ctx.velem = zeros(ctx.nelem,1);
for el = 1:ctx.nelem
    Xe = coords(:, connect(1:4,el));
    ctx.velem(el) = abs(det([Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1), Xe(:,4)-Xe(:,1)]))/6;
end

% -------------------------
% 4) Viscoelastic material
% -------------------------
[ctx.M, ctx.Eo, ctx.Em, ctx.rhom, ctx.tEi] = ViscMod(C.mat);

% -------------------------
% 5) Cohesive law parameters
% -------------------------
ctx.sc1 = C.sc1;
ctx.a1  = C.a1;
ctx.a2  = C.a2;
ctx.Dym = C.Dym;

% -------------------------
% 6) Face operator R and pattern indq0 (matches F0a_core)
% -------------------------
ctx.R = [6 -1 -1 0 -4 0 
        -1 6 -1 0 0 -4 
        -1 -1 6 -4 0 0 
         0 0 -4 32 16 16 
        -4 0 0 16 32 16 
         0 -4 0 16 16 32]/180;

ctx.indq0 = reshape(reshape(1:36,6,6).',36,1);
ctx.nfnodes  = 6;
ctx.nfnodes2 = 36;
ctx.ncoord = 3;

% -------------------------
% 7) Legacy geometric masks + BCs
% -------------------------
e    = 1e-8;
b1   = 1.5;
lcoh = 3;

% continuation of crack plane region (used to RESTRICT cohesive faces)
cind = find(coords(3,:)<e & coords(2,:)<b1+lcoh+e);
state.cind = cind;

% fixed sets (legacy)
fixed1 = find(coords(3,:)<e & coords(2,:)>b1+lcoh-e); % fix z only
fixed2 = find(coords(2,:)<e);                          % fix y only
fixed3 = find(coords(1,:)<e);                          % fix x only

% DOF indices: x=3*i-2, y=3*i-1, z=3*i
ctx.fixvars = sort([3*fixed3(:)-2; 3*fixed2(:)-1; 3*fixed1(:)]);

% crack-front control line (legacy cohsym)
cohsym = find(coords(3,:)<e & coords(1,:)<e & coords(2,:)<b1+lcoh+e).';
tmp = sortrows([cohsym, coords(2,cohsym)'], 2);
cohsym = tmp(:,1);

state.vc = 3*cohsym(:);  % z-DOFs
% ci0 should be based on y-coordinate, not node id
state.ci0 = find(coords(2,cohsym) > b1+lcoh-e, 1, 'first');

% top surface nodes for debug (legacy eind)
d3 = max(coords(3,:));
state.eind = find(abs(coords(3,:) - d3) < e);

% -------------------------
% 8) Boundary faces detection (unique corner-triplet)
% -------------------------
faces = [1 2 3 5 6 7;
         1 4 2 8 9 5;
         2 4 3 9 10 6;
         3 4 1 10 8 7];

allKey = zeros(3,4*ctx.nelem);
allE   = zeros(1,4*ctx.nelem);
allF   = zeros(1,4*ctx.nelem);
k = 0;

for el = 1:ctx.nelem
    conn = connect(:,el);
    for f = 1:4
        k = k+1;
        fn = conn(faces(f,1:3));
        allKey(:,k) = sort(fn(:));
        allE(k) = el;
        allF(k) = f;
    end
end

[~,ord] = sortrows(allKey.', [1 2 3]);
keyS = allKey(:,ord);
eS   = allE(ord);
fS   = allF(ord);

dupPrev = [false, all(keyS(:,2:end)==keyS(:,1:end-1),1)];
dupNext = [all(keyS(:,2:end)==keyS(:,1:end-1),1), false];
isBoundary = ~(dupPrev | dupNext);

bE = eS(isBoundary);
bF = fS(isBoundary);

triArea = @(p1,p2,p3) 0.5*norm(cross(p2-p1, p3-p1));

% For restricting cohesive faces to the intended region:
cindMask = false(1, ctx.nnodes);
cindMask(cind) = true;

cohFaces  = zeros(0,2);
loadFaces = zeros(0,2);
cohArea   = zeros(0,1);
loadArea  = zeros(0,1);

for i = 1:numel(bE)
    el = bE(i); f = bF(i);
    conn = connect(:,el);
    fn6 = conn(faces(f,:));
    cn  = fn6(1:3);
    zc  = coords(3,cn);

    % cohesive faces: z≈0 AND (corner nodes) belong to cind (y<b1+lcoh region)
    if max(abs(zc - 0)) < e && all(cindMask(cn))
        cohFaces(end+1,:) = [el f]; %#ok<AGROW>
        p1 = coords(:,cn(1)); p2 = coords(:,cn(2)); p3 = coords(:,cn(3));
        cohArea(end+1,1) = triArea(p1,p2,p3); %#ok<AGROW>

    % loaded faces: z≈d3
    elseif max(abs(zc - d3)) < e
        loadFaces(end+1,:) = [el f]; %#ok<AGROW>
        p1 = coords(:,cn(1)); p2 = coords(:,cn(2)); p3 = coords(:,cn(3));
        loadArea(end+1,1) = triArea(p1,p2,p3); %#ok<AGROW>
    end
end

ctx.mc = size(cohFaces,1);
ctx.me = size(loadFaces,1);

% -------------------------
% 9) Build jc (cohesive dofs) and cfar (areas)
% -------------------------
ctx.cfar = cohArea;
ctx.jc   = zeros(ctx.nfnodes*ctx.mc,1);
for i = 1:ctx.mc
    el = cohFaces(i,1); f = cohFaces(i,2);
    fn6 = connect(faces(f,:), el);
    ctx.jc((i-1)*ctx.nfnodes+(1:ctx.nfnodes)) = 3*fn6(:); % z-DOFs
end

% -------------------------
% 10) Build external load maps je and re (uniform traction in z)
% External traction (constant) on quadratic face: corner weights = 0, midside = A/3
ctx.je = zeros(ctx.me*3,1);
ctx.re = zeros(ctx.me*3,1);
for i = 1:ctx.me
    el = loadFaces(i,1);
    f  = loadFaces(i,2);

    fn6 = connect(faces(f,:), el);   % [6 nodes]
    mid = fn6(4:6);                  % midside nodes

    ind = (i-1)*3 + (1:3);
    ctx.je(ind) = 3*mid(:);          % z-DOFs
    ctx.re(ind) = loadArea(i)/3;     % consistent weights
end

% -------------------------
% 11) Allocate evolving state
% -------------------------
state.U      = zeros(ctx.ndof,1);
state.dU0    = zeros(ctx.ndof,1);
state.Fpsig  = zeros(ctx.ndof,1);
state.Ftsig  = zeros(ctx.ndof,1);
state.S      = zeros(6,6,ctx.nip3*ctx.nelem,ctx.M);
state.sigzip = zeros(ctx.nelem*ctx.nip3,1);
state.t1     = zeros(C.nt+1,1);

state.ci = 1;
state.kD = NaN;
state.sig = C.sig_target;

end
