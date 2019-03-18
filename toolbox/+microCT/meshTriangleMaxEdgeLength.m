function L = meshTriangleMaxEdgeLength(F, V)
%MESHTRIANGLEMAXEDGELENGTH Compute the maximum edge length for each
%triangle of the given triangulation.
%
%  L = meshTriangleMaxEdgeLength(F, V)
%
%    F - Mx3 triangle matrix containing vertex indices
%    V - Nx3 matrix of vertex positions
%
%    L - Mx1 vector of per-triangle maximum edge length
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.


% circular expansion of F
F = [F, F(:, 1)];

L = zeros(size(F, 1), 1);

for ii = 1:3
  
  L = max(L, sum((V(F(:, ii + 1), :) - V(F(:, ii), :)).^2, 2));
  
end

L = sqrt(L);

end

