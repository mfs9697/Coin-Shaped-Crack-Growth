function [ctx, state] = init_Notch3D02(C)
%INIT_NOTCH3D02  Build ctx/state for Notch3D case (globals-free).
%
% This init loads Geometry/crd05.txt and con05.txt (T4), converts to T10,
% builds boundary faces, cohesive faces and loaded faces maps (jc/je/re/cfar),
% and sets fixvars and vc control DOFs.

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

% Convert to the orientation used in your legacy: coords is 3 x nnodes
coords = coords_in.';                 % 3 x nnodes
connect4 = connect4_in.';             % 4 x nelem

[coords, connect] = T4toT10(coords, connect4);   % connect is 10 x nelem

ctx.coords   = coords;
ctx.connect  = connect;
ctx.nnodes   = size(coords,2);
ctx.nelem    = size(connect,2);
ctx.nelnodes = 10;
ctx.ndof     = 3*ctx.nnodes;
ctx.eldf     = 3*ctx.nelnodes;

% -------------------------
% 2) Integration points / weights
% -------------------------
[ctx.xip3, ctx.w3] = IntegrationParameters();
ctx.nip3 = numel(ctx.w3);

% -------------------------
% 3) Element volumes (needed in some legacy post-proc; harmless to store)
% -------------------------
ctx.velem = zeros(ctx.nelem,1);
for e = 1:ctx.nelem
    Xe = coords(:, connect(1:4,e)); % corners only
    ctx.velem(e) = abs(det([Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1), Xe(:,4)-Xe(:,1)]))/6;
end

% -------------------------
% 4) Viscoelastic material (keep consistent with your kernels)
% -------------------------
% If your repo ViscMod() expects no inputs, keep it that way.
[ctx.M, ctx.Eo, ctx.Em, ctx.rhom, ctx.tEi] = ViscMod();

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
ctx.R = [6 -1 -1 0 -4 0 ...
        -1 6 -1 0 0 -4 ...
        -1 -1 6 -4 0 0 ...
         0 0 -4 32 16 16 ...
        -4 0 0 16 32 16 ...
         0 -4 0 16 16 32]/180;

ctx.indq0 = reshape(reshape(1:36,6,6).',36,1);
ctx.nfnodes  = 6;
ctx.nfnodes2 = 36;

% -------------------------
% 7) Boundary faces detection (unique corner-triplet)
% -------------------------
tol = 1e-8;

faces = [1 2 3 5 6 7;
         1 4 2 8 9 5;
         2 4 3 9 10 6;
         3 4 1 10 8 7];    % face nodes on a T10 tetra

allKey = zeros(3,4*ctx.nelem);
allE   = zeros(1,4*ctx.nelem);
allF   = zeros(1,4*ctx.nelem);

k = 0;
for e = 1:ctx.nelem
    conn = connect(:,e);
    for f = 1:4
        k = k+1;
        fn = conn(faces(f,1:3));
        allKey(:,k) = sort(fn(:));
        allE(k) = e;
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

zmax = max(coords(3,:));
ymin = min(coords(2,:));

cohFaces = zeros(0,2);
loadFaces = zeros(0,2);
cohArea  = zeros(0,1);
loadArea = zeros(0,1);

for i = 1:numel(bE)
    e = bE(i); f = bF(i);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));
    cn  = fn6(1:3);
    zc  = coords(3,cn);

    if max(abs(zc - 0)) < tol
        cohFaces(end+1,:) = [e f]; %#ok<AGROW>
        p1 = coords(:,cn(1)); p2 = coords(:,cn(2)); p3 = coords(:,cn(3));
        cohArea(end+1,1) = triArea(p1,p2,p3); %#ok<AGROW>
    elseif max(abs(zc - zmax)) < tol
        loadFaces(end+1,:) = [e f]; %#ok<AGROW>
        p1 = coords(:,cn(1)); p2 = coords(:,cn(2)); p3 = coords(:,cn(3));
        loadArea(end+1,1) = triArea(p1,p2,p3); %#ok<AGROW>
    end
end

ctx.mc = size(cohFaces,1);
ctx.me = size(loadFaces,1);

% -------------------------
% 8) Build jc (cohesive dofs) and cfar (areas)
% Using z-DOFs as normal opening component.
% -------------------------
ctx.cfar = cohArea;
ctx.jc   = zeros(ctx.nfnodes*ctx.mc,1);
for i = 1:ctx.mc
    e = cohFaces(i,1); f = cohFaces(i,2);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));
    ctx.jc((i-1)*ctx.nfnodes+(1:ctx.nfnodes)) = 3*fn6(:); % z-DOFs
end

% -------------------------
% 9) Build external load maps je and re (uniform traction in z)
% -------------------------
ctx.je = zeros(ctx.nfnodes*ctx.me,1);
ctx.re = zeros(ctx.nfnodes*ctx.me,1);
for i = 1:ctx.me
    e = loadFaces(i,1); f = loadFaces(i,2);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));
    ind = (i-1)*ctx.nfnodes+(1:ctx.nfnodes);
    ctx.je(ind) = 3*fn6(:);               % z-DOFs
    ctx.re(ind) = loadArea(i)/ctx.nfnodes; % simple consistent distribution
end

% -------------------------
% 10) Fixvars: clamp y=ymin nodes (all 3 DOFs fixed)
% -------------------------
fixedNodes = find(abs(coords(2,:) - ymin) < tol);
ctx.fixvars = sort([3*fixedNodes-2, 3*fixedNodes-1, 3*fixedNodes].');

% -------------------------
% 11) state.vc: control DOFs along x≈0, z≈0 sorted by y
% -------------------------
frontNodes = find(abs(coords(1,:))<tol & abs(coords(3,:))<tol).';
[~,ordF] = sort(coords(2,frontNodes));
frontNodes = frontNodes(ordF);
state.vc = 3*frontNodes(:);  % z-DOFs

% -------------------------
% 12) Allocate evolving state
% -------------------------
state.U      = zeros(ctx.ndof,1);
state.dU0    = zeros(ctx.ndof,1);
state.Fpsig  = zeros(ctx.ndof,1);
state.Ftsig  = zeros(ctx.ndof,1);
state.S      = zeros(6,6,ctx.nip3*ctx.nelem,ctx.M);
state.sigzip = zeros(ctx.nelem*ctx.nip3,1);

state.t1     = zeros(C.nt+1,1);

state.ci = 1;
state.kD = NaN;        % set in initial elastic step
state.sig = C.sig_target;  % used as initial guess level in some solves
end
