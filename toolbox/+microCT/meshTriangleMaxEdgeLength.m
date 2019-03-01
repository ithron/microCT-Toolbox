function L = meshTriangleMaxEdgeLength(F, V)
%MESHTRIANGLEMAXEDGELENGTH Summary of this function goes here
%   Detailed explanation goes here


% circular expansion of F
F = [F, F(:, 1)];

L = zeros(size(F, 1), 1);

for ii = 1:3
  
  L = max(L, sum((V(F(:, ii + 1), :) - V(F(:, ii), :)).^2, 2));
  
end

L = sqrt(L);

end

