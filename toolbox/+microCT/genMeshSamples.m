function [Vnew, UVnew, Inew, Nnew, CtThnew] = genMeshSamples(F, V, CtTh)
%GENMESHSAMPLES Summary of this function goes here
%   Detailed explanation goes here

TR = triangulation(F, V);

N = TR.vertexNormal;

Subdivs = floor(microCT.meshTriangleMaxEdgeLength(F, V)) + 1;

UniqueSubdivs = unique(Subdivs);

Vnew = [];
Nnew = [];
UVnew = [];
Inew = [];

for si=1:length(UniqueSubdivs)
  
  sdiv = UniqueSubdivs(si);
  FL = Subdivs == sdiv;
  FI = find(FL);
  
  [Vi, UVi, Ii, Ni] = microCT.subsampleTriangulation(V, F(FL, :), sdiv, N);
  
  Ii = FI(Ii);
  
  Vnew = [Vnew; Vi];
  Nnew = [Nnew; Ni];
  UVnew = [UVnew; UVi];
  Inew = [Inew; Ii];
  
end

[Vnew, IA, ~] = uniquetol(Vnew, 0.5/max(abs(V(:))), 'ByRows', true);
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

