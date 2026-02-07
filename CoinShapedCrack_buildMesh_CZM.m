function model = CoinShapedCrack_buildMesh_CZM( ...
    L, H, a0, a2, ...
    ncoh, hgradOut, ...
    cohBndID, meshFile)
% CoinShapedCrack_buildMesh_CZM_onlyOuter
% -------------------------------------------------------------------------
% COMSOL mesh builder for 1/8 coin-shaped crack block (CZM only),
% with a single refinement annulus on z=0 between radii a0 and a2.
%
% INPUT:
%   L, H       - block dimensions (m)
%   a0         - crack radius (m)
%   a2         - outer ring radius (m)  (cohesive annulus ends here)
%   ncoh       - ring resolution: hOut = (a2-a0)/ncoh
%   hgradOut   - growth factor from cohesive ring to the bulk
%   cohBndID   - boundary ID (geom dim=2) of the quarter-annulus on z=0
%                (between a0 and a2). If unknown, pass [] and read console.
%   meshFile   - output .mat filename (optional)
%
% OUTPUT:
%   model      - COMSOL model object
%
% SIDE EFFECT:
%   Exports mesh to meshFile (default: CoinShapedCrack_mesh_CZM_raw.mat)

import com.comsol.model.*
import com.comsol.model.util.*

if nargin < 8 || isempty(meshFile)
    meshFile = 'CoinShapedCrack_mesh_CZM_raw.mat';
end

%% Start COMSOL server if needed
try
    mphtags('-show');
catch
    mphstart;
end

fprintf('Building CZM-only mesh (outer ring only)\n');
fprintf('Radii: a0 = %.6g, a2 = %.6g\n', a0, a2);

% numeric ring size
hOut = (a2 - a0) / max(ncoh, 1);
fprintf('Outer ring: hOut = %.4e m (ncoh=%g), hgradOut=%g\n', hOut, ncoh, hgradOut);

%% Create model
model = ModelUtil.create('CoinCrack_CZM_onlyOuter');
comp  = model.component.create('comp1', true);
geom1 = comp.geom.create('geom1', 3);
mesh1 = comp.mesh.create('mesh1');

% optional physics (safe to keep; doesn't hurt)
comp.physics.create('solid', 'SolidMechanics', 'geom1');

%% Parameters
model.param.set('L',  sprintf('%g[m]', L));
model.param.set('H',  sprintf('%g[m]', H));
model.param.set('a0', sprintf('%g[m]', a0));
model.param.set('a2', sprintf('%g[m]', a2));

%% Geometry: block 0<=x<=L, 0<=y<=L, 0<=z<=H
blk1 = geom1.create('blk1', 'Block');
blk1.set('size', {'L','L','H'});
blk1.set('base', 'corner');
blk1.set('pos',  {'0','0','0'});

%% Work plane at z=0 with two circles (a0, a2)
wp1 = geom1.create('wp1', 'WorkPlane');
wp1.set('quickplane', 'xy');     % z=0
wp1g = wp1.geom;

% crack front radius a0
c0 = wp1g.create('c0', 'Circle');
c0.set('r',     'a0');
c0.set('angle', 90);

% outer radius a2
c2 = wp1g.create('c2', 'Circle');
c2.set('r',     'a2');
c2.set('angle', 90);

wp1g.run;
geom1.run;
geom1.run('fin');

%% Cohesive annulus boundary selection (explicit by ID)
% IMPORTANT: this ID depends on the geometry build, and will likely differ
% from your previous (INNER_ID, OUTER_ID) pair.
if isempty(cohBndID)
    warning(['cohBndID is empty. You must provide the boundary ID of the ', ...
             'quarter-annulus (a0..a2) on z=0.\n', ...
             'Tip: open the model in COMSOL GUI and inspect boundary IDs, ', ...
             'or temporarily print selections using mphviewselection / GUI.\n', ...
             'For now no local refinement will be applied (global mesh only).']);
else
    comp.selection.create('bCoh', 'Explicit');
    comp.selection('bCoh').geom(2);
    comp.selection('bCoh').set(cohBndID);
    comp.selection('bCoh').label('bCoh_outerRing');
end

%% Mesh: FreeTet with one numeric size feature
ft = mesh1.create('ftet1', 'FreeTet');

% Global/background size
mesh1.feature('size').set('hauto', 9);

% Local size on cohesive annulus
if ~isempty(cohBndID)
    szOut = ft.create('sizeCoh', 'Size');
    szOut.selection.named('bCoh');
    szOut.set('custom',      true);
    szOut.set('hminactive',  true);
    szOut.set('hmaxactive',  true);
    szOut.set('hgradactive', true);
    szOut.set('hmin',  sprintf('%g[m]', hOut));
    szOut.set('hmax',  sprintf('%g[m]', hOut));
    szOut.set('hgrad', sprintf('%g',   hgradOut));
end

mesh1.run;

%% Export mesh to MATLAB
[meshstats, meshdata] = mphmeshstats(model, 'mesh1');

vtx = meshdata.vertex;
X = vtx(1,:).';
Y = vtx(2,:).';
Z = vtx(3,:).';

idxTet = find(strcmp(meshdata.types, 'tet'), 1);
if ~isempty(idxTet)
    Tet0 = double(meshdata.elem{idxTet}.');   % 0-based
    Tet  = Tet0 + 1;                          % 1-based
else
    warning('No tetrahedral elements found.');
    Tet = [];
end

caseType = "CZM";
save(meshFile, 'X','Y','Z','Tet','meshstats','meshdata', ...
               'L','H','a0','a2', ...
               'ncoh','hgradOut','caseType', ...
               'hOut','cohBndID');

fprintf('Mesh exported to %s\n', meshFile);

end
