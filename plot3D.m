function h = plot3D( varargin )
%PDEPLOT3D Plot mesh and solution for 3-D PDE
%   H=PDEPLOT3D(PDEM) plots the mesh for a 3-D PDE defined by the PDEModel
%   object PDEM. The handles to the plotted graphics objects are returned 
%   in the optional output argument H.
%
%   H=PDEPLOT3D(NODES, ELEMENTS) plots the mesh for a 2-D PDE defined by the 
%   NODES and ELEMENTS matrices. NODES is a 3-by-NumNodes matrix representing 
%   nodal coordinates. ELEMENTS is a 4-by-NumElements or 10-by-NumElements 
%   matrix representing element connectivity in terms of node IDs. 
%
%   H=PDEPLOT3D(...,'Name1',Value1, 'Name2',Value2,...) plots the solution
%   to a 3-D PDE over the specified mesh or elements. 
%
%   The Name-Value pairs are used to specify the solution values and
%   plotting styles are as follows:
%
%   Name        Value/{Default}   Description
%   ----------------------------------------------------------------
%   ColormapData            data          - Colormap plot.
%                                           ColormapData is a vector of nodal 
%                                           solution values.
%   Mesh                    {off} | on    - Display mesh edges. This option
%                                           is used in conjunction with
%                                           ColormapData. It is ignored if
%                                           ColormapData it not specified.
%   NodeLabels              {off} | on    - Display node labels on the 
%                                           mesh boundary. Use in
%                                           conjunction with FaceAlpha.
%   ElementLabels           {off} | on    - Display element labels on the 
%                                           mesh boundary. Use in
%                                           conjunction with FaceAlpha.
%   FaceAlpha               data          - Transparency of faces.
%                                           Scalar in the range [0, 1],
%                                           default = 1
%   Deformation             Displacement  - For a structural analysis model, 
%                                           use this option to plot the
%                                           deformed shape. Value must be
%                                           the displacement structure,
%                                           available in
%                                           StaticStructuralResults object.
%  DeformationScaleFactor  data           - Defines the scaling factor for
%                                           plotting deformed shape. Use
%                                           this option in conjunction with
%                                           Deformation to override the
%                                           default. Default value is
%                                           computed internally based on
%                                           the geometry dimension and
%                                           magnitude of deformation.
% 
%
%   Example 1: Plot a 3-D mesh.
%     pdem = createpde();
%     importGeometry(pdem,'Block.stl');
%     generateMesh(pdem,'Hmax',20);
%     pdeplot3D(pdem,'NodeLabels','on', 'FaceAlpha',0.5);
%
%   Example 2: Solve a cantilever beam problem and plot the deformed shape.
%     pdem = createpde('structural','static-solid');
%     importGeometry(pdem,'SquareBeam.STL');
%     pdegplot(pdem,'FaceLabels','on')
%     generateMesh(pdem);
%     structuralProperties(pdem,'PoissonsRatio',0.3,'YoungsModulus',210E3);
%     structuralBC(pdem,'Face',6,'Constraint','fixed');
%     structuralBoundaryLoad(pdem,'Face',5,'SurfaceTraction',[0,0,-2]);
%     R = solve(pdem);
%     % Plot deformed shape with von Mises stress as contours 
%     figure    
%     pdeplot3D(pdem,'ColorMapData',R.VonMisesStress,'Deformation',R.Displacement)
%     % Plot deformed shape with a scale factor of 100 and von Mises Stress as contours
%     figure
%     pdeplot3D(pdem,'ColorMapData',R.VonMisesStress,'Deformation',...
%                   R.Displacement,'DeformationScaleFactor',100)
%
%   See also: CREATEPDE, pde.PDEMODEL, PDECONT, PDEGPLOT, pde.FEMesh, PDESURF

%   Copyright 2014-2018 The MathWorks, Inc.


if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

nargs = nargin;
if nargs < 1
  error(message('pde:pdeplot:nargin'))
end

thepde = [];
nodeAndElementDef = false;
if isa(varargin{1}, 'struct') 
    themesh=varargin{1};
    p=themesh.Points;
    t=themesh.Elements;
    varargin = {varargin{2:end}};
elseif isa(varargin{1}, 'pde.FEMesh') 
  themsh = varargin{1};
  if rem(nargs,2)~=1
      error(message('pde:pdeplot:NoParamPairs'))
  end
  if size(themsh.Nodes,1) ~= 3
      error(message('pde:pdeplot:NotThreeD'))
  else      
      [p,~,t] = themsh.meshToPet();    
      varargin = {varargin{2:end}};      
  end  
else
    if nargs < 2
        error(message('pde:pdeplot:InvalidArgs'))
    elseif rem(nargs,2)~=0
        error(message('pde:pdeplot:NoParamPairs'))
    end
    p = varargin{1};
    t = varargin{2};
    nodeAndElementDef = true;
    validateattributes(varargin{1},{'numeric'},{'real', 'finite', 'nonsparse', 'nonnan'},'pdeplot3D',inputname(1));
    validateattributes(varargin{2},{'numeric'},{'real', 'integer'},'pdeplot3D',inputname(2));
    varargin = {varargin{3:end}};    
end

if size(p,1) ~= 3
   error(message('pde:pdeplot:NotThreeD'))
end   

parser = inputParser;
addParameter(parser,'colormapdata', [], @isValidU);
addParameter(parser,'flowdata', [], @isnumeric);
addParameter(parser,'Deformation', []);
addParameter(parser,'DeformationScaleFactor', [], @isnumeric);
addParameter(parser,'NodeLabels', 'off', @isValidNdLabelOption);
addParameter(parser,'ElementLabels', 'off', @isValidElLabelOption);
addParameter(parser,'Mesh', 'off', @isValidMeshDispOption);
addParameter(parser,'FaceAlpha', 1, @isnumeric);
addParameter(parser,'EdgeColor', 'k');
addParameter(parser,'FaceColor', 'cyan');
parse(parser,varargin{:});
colormapdata = parser.Results.colormapdata;
flowdata = parser.Results.flowdata;
faceAlpha = parser.Results.FaceAlpha;
showMesh = parser.Results.Mesh;
showNodeLabels = parser.Results.NodeLabels;
showElemLabels = parser.Results.ElementLabels;
deformation = parser.Results.Deformation;
scaleFactor = parser.Results.DeformationScaleFactor;
faceColor = parser.Results.FaceColor;
edgeColor = parser.Results.EdgeColor;

if(strcmpi(showElemLabels, 'on') && nodeAndElementDef && size(p,1) == 3) 
    error(message('pde:pdeplot:ElemLabelsNeedsCompleteMesh'));
end

if(size(t,1) == 5 || size(t,1) == 11)
    t(end,:)=[];
end   

numElemNodes = size(t,1);
if(numElemNodes ~= 4 && numElemNodes ~= 10 && numElemNodes ~= 20)
  error(message('pde:pdeModel:invalidT', numElemNodes));
end
numNodes = size(p,2);
if ~isempty(colormapdata)
  if(length(colormapdata) ~= numNodes)
    error(message('pde:pdeplot:colormapdataLength'));
  end
end

ha = newplot;
hf = get(ha,'Parent');
set(hf,'Color','white');
set(ha,'ClippingStyle','rectangle');
hold on;

bbox = [min(p(1,:)) max(p(1,:));
        min(p(2,:)) max(p(2,:));
        min(p(3,:)) max(p(3,:))];
    
if ~isempty(deformation) && isa(thepde,'pde.StructuralModel')
    if isempty(scaleFactor)
        MaxDeformationMag = max(sqrt(deformation.ux.^2+deformation.uy.^2+deformation.uz.^2));
        if(MaxDeformationMag > eps())
            scaleFactor = min(bbox(:,2) - bbox(:,1)) /MaxDeformationMag; % Based on lowest bounding box dimension
        else
            scaleFactor = 0;
        end
    end
    p(1,:) = p(1,:) +  scaleFactor*deformation.ux';
    p(2,:) = p(2,:) +  scaleFactor*deformation.uy';
    p(3,:) = p(3,:) +  scaleFactor*deformation.uz';
end

bbox = [min(p(1,:)) max(p(1,:));
        min(p(2,:)) max(p(2,:));
        min(p(3,:)) max(p(3,:))];
 
plotAxisTriad(bbox);


if ~isempty(colormapdata)
    if numElemNodes == 10
        ntri=tetBoundaryFacets(p,t);
        ltri=splitQuadraticTri(ntri);
        j2=size(ntri,1);
    elseif numElemNodes == 20
        face=[1,2,3,5,6,7,8,9,10,17
          4,2,1,14,13,6,5,11,12,18
          4,3,2,16,15,8,7,13,14,19
          4,1,3,12,11,10,9,15,16,20];
        C=unique(t(17:20,:)); j2=0; ntri=zeros(2000,10);
        for i=1:length(C)
            [row,col]=find(t(17:20,:)==C(i));
            if length(row)==1
                j2=j2+1;
                ntri(j2,:)=t(face(row,:),col);
            end
        end
        ntri=ntri(1:j2,:);
        ltri=splitCubicTri(ntri);
    end
    h1=colorbar; ht = hgtransform; colormap('jet');
    if min(colormapdata) ~= max(colormapdata)
        caxis([min(colormapdata) max(colormapdata)]);
    end
  
    h2=patch('Faces',ltri, 'Vertices', p',...
        'FaceVertexCData', colormapdata(:), ...
        'AmbientStrength', .75,  ...
        'EdgeColor', 'none','FaceColor', 'interp', 'parent',ht,...
        'FaceAlpha',faceAlpha,'Clipping','off');
    
  if strcmpi(showMesh, 'on')   
     for i=1:j2
         if numElemNodes==10, ind=ntri(i,[1,4,2,5,3,6,1]);
         else, ind=ntri(i,[1,4,5,2,6,7,3,8,9,1]);
         end
         plot3(p(1,ind),...
               p(2,ind),...
               p(3,ind),'Color',.6*[1,1,1],'Clipping','off');
     end
    %{
    [ltri, lp] = tetBoundaryFacets(p,t(1:4,:));  
    tr = triangulation(ltri, lp);
    tre = (tr.edges())';
    x = lp(:,1);
    y = lp(:,2);
    z = lp(:,3);    
    plot3(x(tre),y(tre),z(tre),'-k');
    %}
      
  end
  if(strcmpi(showNodeLabels, 'on'))
    plotNodeLabels(p, t);
  end
  if(strcmpi(showElemLabels, 'on') && ~isempty(thepde))      
      plotElementLabels(thepde.Mesh);
  end     

  if nargout==1
    h = [h1 h2];
  end
elseif ~isempty(flowdata)
    if(length(flowdata(:)) ~= 3*size(p,2))
       error('Length of flowdata must be 3*number of points');
    end
    if(isvector(flowdata)==1)
      flowdata = reshape(flowdata, numNodes, 3);
    end
    quiver3(p(1,:)', p(2,:)', p(3,:)', flowdata(:,1), flowdata(:,2), flowdata(:,3));
else
  %Plot the boundary triangulation:
  [ltri, lp] = tetBoundaryFacets(p,t(1:4,:));  
  set(gcf, 'renderer', 'opengl');
  hh=trisurf(ltri, lp(:,1),lp(:,2),lp(:,3), ...
    'FaceColor',faceColor, 'EdgeColor',edgeColor,'FaceAlpha', faceAlpha, 'Clipping','off');
  if nargout==1
    h = hh;
  end    
  if(strcmpi(showNodeLabels, 'on'))
    plotNodeLabels(p,t);
  end
  if(strcmpi(showElemLabels, 'on'))
     plotElementLabels(themsh);
  end     
end
axis equal; axis tight; axis off; view(130,20);
hold off;
end


function plotNodeLabels(p,t)
    tri = tetBoundaryFacets(p,t);      
    vispts = unique(tri(:));   
    vislabels = num2str(vispts,'n%d'); 
    warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');    
    tr = triangulation(tri(:,1:3),p');
    warning(warnState);
    [~, ric] = tr.incenter();  
    ric = mean(ric);
    vn = tr.vertexNormal(vispts);
    vn = (vn.*(0.15*ric))';
    px = p(1,vispts)+ vn(1,:);
    py = p(2,vispts)+ vn(2,:);
    pz = p(3,vispts)+ vn(3,:);
    text(px,py,pz, vislabels,'FontSize',8, 'HorizontalAlignment','center', 'Clipping','on');
end



function plotElementLabels(msh)
    ma = msh.MeshAssociation;
    fa = ma.FaceAssociativity;
    facets = cell2mat(fa');
    flabels = num2str(facets(1,:)','e%d'); 
    n = msh.Nodes';
    e = msh.Elements(1:4,:)';
    lfc = ma.FaceCodes(:,1:3); % Linear Face codes
    facets = facets';
    i = sub2ind(size(e), facets(:,1), lfc(facets(:,2),1));
    j = sub2ind(size(e), facets(:,1), lfc(facets(:,2),2));
    k = sub2ind(size(e), facets(:,1), lfc(facets(:,2),3));
    ei = e(i);
    ej = e(j);
    ek = e(k);
    nf = [ei(:), ej(:), ek(:)];
    warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');    
    tr = triangulation(nf,n);
    warning(warnState);
    [ic, ric] = tr.incenter();
    fn = faceNormal(tr);
    ic = ic + (fn.*ric)*0.15;
    text(ic(:,1),ic(:,2),ic(:,3), flabels,'FontSize',8, 'HorizontalAlignment','center', 'Clipping','on');
end

function ok = isValidU(udata)
    if ~isreal(udata) || issparse(udata) || any(~isfinite(udata))
        error(message('pde:pdeplot:InvalidColormapdata'));   
    end     
    ok = true;
end

function ok = isValidMeshDispOption(dispopt)
    if ~ischar(dispopt)
       error(message('pde:pdeplot:GenericMustBeOffOn','mesh'));   
    end
    validatestring(dispopt,{'on', 'off'}) ;
    ok = true;
end

function ok = isValidNdLabelOption(labelopt)
    if ~ischar(labelopt)
       error(message('pde:pdeplot:GenericMustBeOffOn', 'NodeLabels'));   
    end
    validatestring(labelopt,{'on', 'off'}) ;
    ok = true;
end

function ok = isValidElLabelOption(labelopt)
    if ~ischar(labelopt)
       error(message('pde:pdeplot:GenericMustBeOffOn', 'ElementLabels'));   
    end
    validatestring(labelopt,{'on', 'off'}) ;
    ok = true;
end

function t4=splitQuadraticTri(t)
t4Nodes = [1 4 6; 4 5 6; 4 2 5; 6 5 3];
t4 = [t(:,t4Nodes(1,:)); t(:,t4Nodes(2,:)); t(:,t4Nodes(3,:)); t(:,t4Nodes(4,:))];
end

function plotAxisTriad(bbox)
% PLOTAXISTRIAD - Plots cartesian axis triad based on input bounding box
%

%   Copyright 2014-2016 The MathWorks, Inc.

maxdim = max([abs(bbox(1,1)-bbox(1,2)), abs(bbox(2,1)-bbox(2,2)), abs(bbox(3,1)-bbox(3,2))]);
orgn = (ones(3,3).*(bbox(:,1) - 0.3*maxdim))';
vMag = .3*maxdim;
vec = vMag*eye(3);
quiver3(orgn(:,1),orgn(:,2),orgn(:,3),vec(:,1),vec(:,2),vec(:,3), 'Color', 'red');
labPts = (orgn + vec);
text(labPts(:,1), labPts(:,2), labPts(:,3), ['x' 'y' 'z']');
end

