function [coords, connect, edget] = T4to10(coords4, connect4)
%T4TO10  Convert linear T4 tetra mesh to quadratic T10 mesh (mid-edge nodes).
%
% INPUT:
%   coords4  : 3 x N4
%   connect4 : 4 x nelem
%
% OUTPUT:
%   coords   : 3 x N10
%   connect  : 10 x nelem   (T10 connectivity)
%   edget    : ne x 3       [n1, nmid, n2] for each unique edge

nelnodes = 10;
nelem    = size(connect4,2);

edge    = [1,2; 2,3; 3,1; 1,4; 2,4; 3,4];
midedge = 5:10;

connect = zeros(nelnodes, nelem);
connect(1:4,:) = connect4;

% collect all edges (6 per element)
edges = zeros(6*nelem, 2);
edgessort = zeros(6*nelem, 2);

for e = 1:nelem
    for k = 1:6
        ind = (e-1)*6 + k;
        edges(ind,:) = connect4(edge(k,:), e);
        edgessort(ind,:) = sort(edges(ind,:));
    end
end

% unique edges (by sorted node ids)
[~,~,iedge] = unique(edgessort, 'rows');
ne = max(iedge);

n0 = max(connect4, [], 'all');          % number of original nodes

coords = zeros(3, n0 + ne);
coords(:,1:n0) = coords4(:,1:n0);

edget = zeros(ne, 3);

% create mid-edge nodes
n = n0;
for i = 1:ne
    edgei = find(iedge == i);          % occurrences of this unique edge

    % decode which element and which local edge number each occurrence is
    edgeilmn = floor((edgei-1)/6) + 1;
    edgeinum = mod(edgei,6);
    edgeinum(edgeinum==0) = 6;

    n = n + 1;                         % new mid-node index

    % assign mid-node into each element that contains this edge
    for m = 1:numel(edgei)
        connect(midedge(edgeinum(m)), edgeilmn(m)) = n;
    end

    % store [n1, nmid, n2] once
    edget(i,:) = [edges(edgei(1),1), n, edges(edgei(1),2)];

    % place mid-node coordinates
    a = coords(:, edges(edgei(1),1));
    b = coords(:, edges(edgei(1),2));
    coords(:, n) = 0.5*(a + b);
end
end
