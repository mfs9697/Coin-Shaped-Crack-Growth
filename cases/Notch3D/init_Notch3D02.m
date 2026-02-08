function [ctx, state] = init_Notch3D02(cfg)
%INIT_NOTCH3D02  Initialize Notch3D ctx/state (globals-free path).

ctx = struct();
state = struct();

% --- geometry paths ---
here = fileparts(mfilename('fullpath'));
geomDir = fullfile(here,'Geometry');

crdFile = fullfile(geomDir,'crd05.txt');
conFile = fullfile(geomDir,'con05.txt');

coords   = readmatrix(crdFile);
connect4 = readmatrix(conFile);

% --- convert T4 -> T10 ---
[coords, connect] = T4toT10(coords, connect4);

ctx.coords    = coords;
ctx.connect   = connect;

ctx.nnodes    = size(coords,1);
ctx.nelem     = size(connect,2);
ctx.nelnodes  = 10;
ctx.ndof      = 3*ctx.nnodes;
ctx.eldf      = 3*ctx.nelnodes;

% --- integration ---
[ctx.xip3, ctx.w3] = IntegrationParameters();
ctx.nip3 = numel(ctx.w3);

% --- element volumes ---
ctx.velem = zeros(ctx.nelem,1);
for e = 1:ctx.nelem
    Xe = coords(connect(:,e),:);
    ctx.velem(e) = abs(det([Xe(2,:)-Xe(1,:); Xe(3,:)-Xe(1,:); Xe(4,:)-Xe(1,:)]))/6;
end

state.U      = zeros(ctx.ndof,1);
state.dU0    = zeros(ctx.ndof,1);
state.Fpsig  = zeros(ctx.ndof,1);
state.Ftsig  = zeros(ctx.ndof,1);
state.S      = zeros(6,6,ctx.nip3*ctx.nelem,ctx.M);
state.sigzip = zeros(ctx.nelem*ctx.nip3,1);

% ============================================================
% Boundary faces: detect unique faces, then classify:
%   - cohesive faces on z=0
%   - loaded faces on z=zmax
% ============================================================

tol = 1e-8;  % geometric tolerance (tune if needed)

coords = ctx.coords;    % [nnodes x 3] OR [3 x nnodes]? ensure consistent below
connect = ctx.connect;  % [10 x nelem] in your project after T4->T10

% Ensure coords is [3 x nnodes] like your legacy uses
if size(coords,1) ~= 3
    coords = coords.';  % now 3 x nnodes
end
ctx.coords = coords;

nnodes = ctx.nnodes;
nelem  = ctx.nelem;

% Face node numbering for T10 tetra faces (matches your legacy)
faces = [1 2 3 5 6 7;
         1 4 2 8 9 5;
         2 4 3 9 10 6;
         3 4 1 10 8 7];         % 4 faces, 6 nodes each

nfnodes  = 6;
nfnodes2 = 36;
ctx.nfnodes  = nfnodes;
ctx.nfnodes2 = nfnodes2;

% Face operator R (from legacy)
ctx.R = [6 -1 -1 0 -4 0 ...
        -1 6 -1 0 0 -4 ...
        -1 -1 6 -4 0 0 ...
         0 0 -4 32 16 16 ...
        -4 0 0 16 32 16 ...
         0 -4 0 16 16 32]/180;

ctx.indq0 = reshape(reshape(1:36,6,6).',36,1);

% ---- build list of all element faces (by corner nodes as key) ----
% We key faces by their 3 corner nodes (nodes 1:4 in T10 are corners).
allKey = zeros(3, 4*nelem);      % sorted corner ids
allE   = zeros(1, 4*nelem);      % element id
allF   = zeros(1, 4*nelem);      % local face id 1..4
k = 0;
for e = 1:nelem
    conn = connect(:,e);
    for f = 1:4
        fn = conn(faces(f,1:3));         % corner nodes of this face
        k = k+1;
        allKey(:,k) = sort(fn(:));
        allE(k) = e;
        allF(k) = f;
    end
end

% ---- find boundary faces: those with unique corner triplet ----
% Sort keys lexicographically
[~,ord] = sortrows(allKey.', [1 2 3]);
allKeyS = allKey(:,ord);
allES   = allE(ord);
allFS   = allF(ord);

% Count duplicates
isDupPrev = [false, all(allKeyS(:,2:end)==allKeyS(:,1:end-1),1)];
isDupNext = [all(allKeyS(:,2:end)==allKeyS(:,1:end-1),1), false];

isBoundary = ~(isDupPrev | isDupNext);   % appears only once
bE = allES(isBoundary);
bF = allFS(isBoundary);

% ---- helper: triangle area from 3 corner nodes ----
triArea = @(p1,p2,p3) 0.5*norm(cross(p2-p1, p3-p1));

% ---- classify boundary faces ----
zmax = max(coords(3,:));
zmin = min(coords(3,:)); %#ok<NASGU>
ymin = min(coords(2,:));

cohFaces = [];     % rows: [e f]
loadFaces = [];    % rows: [e f]
cohArea  = [];
loadArea = [];

for i = 1:numel(bE)
    e = bE(i); f = bF(i);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));      % 6 nodes on this face

    % use corner nodes to test plane
    cn = fn6(1:3);
    zc = coords(3,cn);

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

ctx.mc = size(cohFaces,1);   % number of cohesive faces
ctx.me = size(loadFaces,1);  % number of loaded faces

% ============================================================
% Build jc (cohesive dof indices) and cfar (areas)
% In your formulation, jc indexes the NORMAL opening component.
% We assume normal is z-direction -> use DOF = 3*node.
% ============================================================
ctx.cfar = cohArea;                        % one area per cohesive face
ctx.jc   = zeros(nfnodes*ctx.mc,1);

for i = 1:ctx.mc
    e = cohFaces(i,1); f = cohFaces(i,2);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));
    ctx.jc((i-1)*nfnodes+(1:nfnodes)) = 3*fn6(:);   % z-DOFs
end

% ============================================================
% Build external load mapping je and re
% Apply uniform traction in z-direction on z=zmax faces.
% Simple consistent nodal load: area/6 per face node (quadratic triangle).
% ============================================================
ctx.je = zeros(nfnodes*ctx.me,1);
ctx.re = zeros(nfnodes*ctx.me,1);

for i = 1:ctx.me
    e = loadFaces(i,1); f = loadFaces(i,2);
    conn = connect(:,e);
    fn6 = conn(faces(f,:));

    ind = (i-1)*nfnodes+(1:nfnodes);
    ctx.je(ind) = 3*fn6(:);              % z-DOFs
    ctx.re(ind) = loadArea(i)/nfnodes;   % distribute uniformly
end

% ============================================================
% Fixvars (Dirichlet DOFs) — clamp nodes at y=ymin
% ============================================================
fixedNodes = find(abs(coords(2,:) - ymin) < tol);
ctx.fixvars = sort([3*fixedNodes-2, 3*fixedNodes-1, 3*fixedNodes].');  % x,y,z

% ============================================================
% state.vc: control DOF indices along crack-front line
% Legacy used abs(x)<e & abs(z)<e and sorted by y.
% We keep z-DOF control: vc(ci) is a single dof index per front position.
% ============================================================
frontNodes = find(abs(coords(1,:))<tol & abs(coords(3,:))<tol).';
[~,ordF] = sort(coords(2,frontNodes));
frontNodes = frontNodes(ordF);

state.vc = 3*frontNodes(:);  % z-DOFs as control coordinates

% Provide defaults for control indices
state.ci = 1;
if ~isfield(state,'kD') || ~isfinite(state.kD)
    state.kD = 0;  % will be set after initial elastic solve
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
