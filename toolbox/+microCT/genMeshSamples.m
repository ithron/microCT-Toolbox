function [Vnew, UVnew, Inew, Nnew, CtThnew] = genMeshSamples(F, V, CtTh)
%GENMESHSAMPLES Densly subsample the given mesh so that there is at least
%one sample per voxel intersecting the mesh.
%
%  [Vnew, UVnew, Inew, Nnew, CtThnew] = genMeshSamples(F, V, CtTh)
%
%    F - Mx3 per triangle vertex index matrix
%    V - Nx3 vertex position matrix
%    CtTh - Nx1 vector of per-vertex cortical thickness values
%
%    Vnew - Kx3 matrix of sample positions
%    UVnew - Kx2 matrix of barycentric sample positions
%    Inew - Kx1 matrix of per-sample triangle index
%    Nnew - Kx3 matrix of per-sample normals
%    CtThnew - Kx1 vector of per-sample cortical thickness values
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

TR = triangulation(double(F), double(V));

N = TR.vertexNormal;

Subdivs = floor(microCT.meshTriangleMaxEdgeLength(F, V)) + 1;

UniqueSubdivs = unique(Subdivs);

nSubdivs = length(UniqueSubdivs);

Vnew = cell(nSubdivs, 1);
Nnew = cell(nSubdivs, 1);
UVnew = cell(nSubdivs, 1);
Inew = cell(nSubdivs, 1);

parfor si=1:nSubdivs
  
  sdiv = UniqueSubdivs(si);
  FL = Subdivs == sdiv;
  FI = find(FL);
  
  [Vi, UVi, Ii, Ni] = microCT.subsampleTriangulation(V, F(FL, :), sdiv, N);
  
  Ii = FI(Ii);
  
  Vnew{si} = Vi;
  Nnew{si} = Ni;
  UVnew{si} = UVi;
  Inew{si} = Ii;
  
end

Vnew = cell2mat(Vnew);
Nnew = cell2mat(Nnew);
UVnew = cell2mat(UVnew);
Inew = cell2mat(Inew);

% [Vnew, IA, ~] = uniquetol(Vnew, 0.5/max(abs(V(:))), 'ByRows', true);
[Vnew, IA] = microCT.filterClosePointsQuick(Vnew, 0.5);
Nnew = Nnew(IA, :);
UVnew = UVnew(IA, :);
Inew = Inew(IA);

% Normalize normals
Nnew = Nnew ./ sqrt(sum(Nnew.^2, 2));

% Compute per sample CtTh
CtThnew = ...
  CtTh(F(Inew, 1)) .* UVnew(:, 1) + ...
  CtTh(F(Inew, 2)) .* UVnew(:, 2) + ...
  CtTh(F(Inew, 3)) .* (1 - sum(UVnew, 2));

end

