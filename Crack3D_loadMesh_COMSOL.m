 function [coords, connect4, model] = Crack3D_loadMesh_COMSOL(meshFile)
% Crack3D_loadMesh_COMSOL
% -------------------------------------------------------------------------
% Mimics the mesh generation block of the original code, but loads 
% the mesh from a COMSOL-generated .mat file instead of building it 
% with addVertex/generateMesh.
%
% INPUT:
%   meshFile - Filename (e.g., 'CoinShapedCrack_mesh_CZM_raw.mat')
%
% OUTPUT:
%   coords   - 3xN array of node coordinates (compatible with solver)
%   connect4 - 4xM array of element connectivity (compatible with solver)
%   model    - The MATLAB PDE model object (for plotting/validation)
% -------------------------------------------------------------------------

fprintf('--- Loading COMSOL Mesh into MATLAB PDE Model ---\n');

% 1. Load the COMSOL Mesh Data
if nargin < 1
    meshFile = 'CoinShapedCrack_mesh_CZM_raw.mat';
end

if ~isfile(meshFile)
    error('Mesh file %s not found. Run the COMSOL builder first.', meshFile);
end

data = load(meshFile);
X = data.X;
Y = data.Y;
Z = data.Z;
Tet = data.Tet; % Mx4 array (1-based indices)

% 2. Create the MATLAB PDE Model Container
model = createpde();

% 3. Import the Mesh into the Model
% geometryFromMesh creates the geometry strictly from the finite elements.
% Nodes must be Transposed to (3, NumNodes) for geometryFromMesh
nodes = [X, Y, Z]'; 
elements = Tet';    % Transpose to (4, NumElements)

geometryFromMesh(model, nodes, elements);

% 4. Mimic the Output of the Old Code
% The old code extracted: coords = model.Mesh.Nodes
% and connect4 = model.Mesh.Elements
coords = model.Mesh.Nodes;
connect4 = model.Mesh.Elements;
connect4 = connect4([1,3,2,4],:);

% Note: The COMSOL mesh is likely already in the positive quadrant [0, L].
% The old code shifted coordinates: coords(1,:) = coords(1,:) + d(1)/2.
% We check if a shift is needed.
min_x = min(coords(1,:));
if min_x < -1e-6
    fprintf('Shifting mesh to positive quadrant...\n');
    coords(1,:) = coords(1,:) + abs(min_x);
    coords(2,:) = coords(2,:) + abs(min(coords(2,:)));
end

% 5. Visualization (Mimics pdeplot3D(model))
figure(1); clf;
pdeplot3D(model);
title(['Imported COMSOL Mesh: ', meshFile], 'Interpreter', 'none');
axis equal;

fprintf('Mesh Loaded: %d Nodes, %d Elements.\n', size(coords,2), size(connect4,2));

    end