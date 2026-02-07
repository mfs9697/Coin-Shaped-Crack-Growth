function [ctx, state] = init_problem(cfg)
%INIT_PROBLEM  Build ctx/state from cfg, load mesh, prepare operators.

% --- mesh / geometry ---
L = cfg.geom.L; H = cfg.geom.H; Rc = cfg.geom.Rc;

if ~isfile(cfg.meta.meshFile)
    % first run: pass [] to force a warning; then set cohBndID after you inspect
    cohBndID = 5; % <-- fill after first GUI inspection
    CoinShapedCrack_buildMesh_CZM(cfg.geom.L, cfg.geom.H, cfg.geom.Rc, cfg.geom.cr_a2, ...
        cfg.mesh.ncoh, cfg.mesh.hgrad_out, cohBndID, cfg.meta.meshFile);
end

[coords, connect4, model] = Crack3D_loadMesh_COMSOL(cfg.meta.meshFile);
nelem = size(connect4,2);

% --- element choice locked ---
nelnodes = cfg.elem.nelnodes;
if nelnodes ~= 10
    error('Locked to quadratic tetrahedra (nelnodes=10).');
end

% face numbering for T10
face = [1,2,3,5,6,7
        4,2,1,9,5,8
        4,3,2,10,6,9
        4,1,3,8,7,10];

meNode = max(connect4,[],[1,2]);
[coords,connect,~] = T4to10(coords(:,1:meNode), connect4);

assert(size(connect,1)==10, 'T4to10 produced %dx%d connectivity (expected 10xnelem).', size(connect,1), size(connect,2));


nnodes = size(coords,2);
ncoord = 3;
nfnodes = 6; nfnodes2 = nfnodes^2;

eldf = ncoord*nelnodes;
ndof = ncoord*nnodes;

nip3 = cfg.elem.nip3;
[xip3,w3,Nextr] = IntegtationParameters(nip3, nelnodes, @N_3DT);

% element volumes
velem = zeros(nelem,1);
for i=1:nelem
    coordc = coords(:,connect(1:4,i));
    velem(i) = det([coordc;1 1 1 1]')/6;
end

% --- cohesive geometry faces (cind), symmetry, fixed dofs ---
e = 1e-8;
tol_radius = 1e-2;
cind = find(abs(coords(3,:)) < e & hypot(coords(1,:), coords(2,:)) > (Rc - tol_radius));

fixed2 = find(abs(coords(1,:))<e);
fixed3 = find(abs(coords(2,:))<e);
fixvars = [ncoord*fixed2'-2; ncoord*fixed3'-1];

% --- load faces (top z=H) ---
d = [L, L, H];
eind = find(abs(coords(3,:)-d(3))<e);

% T10 face operators
R = [ 6 -1 -1  0 -4  0
     -1  6 -1  0  0 -4
     -1 -1  6 -4  0  0
      0  0 -4 32 16 16
     -4  0  0 16 32 16
      0 -4  0 16 16 32]/180;
R0 = [0;0;0;1;1;1]/3;

% find je/re (loaded faces) and jc/cfar (cohesive faces)
je = zeros(1e3*nfnodes,1); re = je; me = 0; efar = zeros(1e3,1);
jc = zeros(1e3*nfnodes,1); mc = 0; cfar = zeros(1e3,1);
crf = zeros(nfnodes,1e3);

for lmn=1:nelem
    ecom = intersect(connect(:,lmn), eind);
    if length(ecom)==nfnodes
        me = me+1;
        for i=1:4
            facenodes = connect(face(i,:), lmn);
            if length(intersect(ecom, facenodes))==nfnodes
                ind = (me-1)*nfnodes + (1:nfnodes);
                je(ind) = ncoord*facenodes;
                efar(me)= polyarea(coords(1,facenodes(1:3)), coords(2,facenodes(1:3)));
                re(ind) = efar(me)*R0;
                break;
            end
        end
    end

    ccom = intersect(connect(:,lmn), cind);
    if length(ccom)==nfnodes
        mc = mc+1;
        for i=1:4
            facenodes = connect(face(i,:), lmn);
            if length(intersect(ccom, facenodes))==nfnodes
                crf(:,mc)= facenodes(1:nfnodes);
                ind = (mc-1)*nfnodes + (1:nfnodes);
                jc(ind) = ncoord*facenodes;
                cfar(mc)= polyarea(coords(1,facenodes(1:3)), coords(2,facenodes(1:3)));
                break;
            end
        end
    end
end

je = je(1:me*nfnodes); re = re(1:me*nfnodes);
jc = jc(1:mc*nfnodes); cfar = cfar(1:mc);
crf = crf([1,4,2,5,3,6],1:mc);

% --- CZM parameters ---
sc1 = cfg.czm.sig_max;
a1  = cfg.czm.tsl_a1;
a2  = cfg.czm.tsl_a2;
c   = (3 - 2*a1 + 3*a2)/6;
Dym = cfg.czm.phi / c / sc1 * 1e-7;

% --- viscoelastic moduli ---
[M,Eo,Em,rhom,tEi] = ViscMod(cfg.mat);

% --- ctx (immutable) ---
ctx = struct();
ctx.coords   = coords; 
ctx.connect  = connect; 
ctx.nelem    = nelem;
ctx.ncoord   = ncoord; 
ctx.nelnodes = nelnodes; 
ctx.ndof     = ndof; 
ctx.eldf     = eldf;
ctx.nfnodes  = nfnodes; 
ctx.nfnodes2 = nfnodes2;
ctx.nip3     = nip3; 
ctx.xip3     = xip3; 
ctx.w3       = w3; 
ctx.Nextr    = Nextr;
ctx.face     = face; 
ctx.velem    = velem;
ctx.fixvars  = fixvars;
ctx.je       = je;
ctx.re       = re; 
ctx.me       = me;
ctx.jc       = jc; 
ctx.cfar     = cfar;
ctx.mc       = mc;
ctx.R        = R;
ctx.indq0    = reshape(reshape(1:nfnodes2,nfnodes,nfnodes)',nfnodes2,1);
ctx.Dym      = Dym; 
ctx.sc1      = sc1; 
ctx.a1       = a1; 
ctx.a2       = a2;
ctx.M        = M;
ctx.Em       = Em;
ctx.rhom     = rhom;
ctx.tEi      = tEi;
ctx.Eo       = Eo;
ctx.E30      = cfg.mat.Et_inst;
ctx.crf      = crf; 
ctx.eind     = eind; 
ctx.d                =d;

% --- state (mutable) ---
state = struct();
state.U = zeros(ndof,1);
state.dU = zeros(ndof,1);
state.dU0 = zeros(ndof,1);
state.Fpsig = zeros(ndof,1);
state.Ftsig = zeros(ndof,1);

state.S = zeros(6,6,nip3*nelem,M);
state.zetamn = rhom;     % initial shape same as your current approach
state.etamn  = rhom;

% crack-front tracking
yb = find(abs(coords(1,:))<e & abs(coords(3,:))<e)';
yb = sortrows([yb,coords(2,yb)'],2);
ci0 = find(yb(:,2)>Rc-1e-8,1,'first');
state.yb = yb;
state.ci = ci0;
state.ci0 = ci0;
state.vc = yb(:,1)*ncoord;
state.kD = 0;

% solver scalars
state.sig = cfg.solve.sig;

% precompute initial stiffness for first step later
state.K = [];
state.tE = [];
end
