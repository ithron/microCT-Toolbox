function Y = accumVertexData(F, UV, I, X)
%ACCUMVERTEXDATA Accumulates densly sampled surface data and compute a
%per-vertex mean value.
%
%  Y = accumVertexData(F, UV, I, X)
%
%    F - Mx3 triangle matrix containing vertex indices
%    UV - Nx2 matrix of sample barycentric coordinates
%    I - Nx1 vector of per-sample triangle indices
%    X - NxK per-sample input matrix
%
%    Y - max(F(:))xK output matrix
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

% Compute weights
W = [UV(:, 1), UV(:, 2), 1 - sum(UV, 2)];

% Copmute weighted inputs
% XW has dimension size(X, 1) x 3 x size(X, 2)
XW = W .* permute(X, [1, 3, 2]);
% reshape to 3 * size(X, 1) x size(X, 2)
XW = reshape(XW, 3 * size(X, 1), size(X, 2));

% For each point get vertex -> point correspondences.
% Since each points lies in one triangle, each input point corresponds to
% 3 vertices, therefore generate 3 index lists pts -> vertex index
IV = [F(I, 1); F(I, 2); F(I, 3)];

% Accumulate all weighted input that correspond to the same vertex
N = max(F(:));
Y = zeros(N, size(X, 2), class(X));

% For each row of the output of accumarray find out with actual vertex
% index it corresponds to
IVperm = accumarray(IV, IV, [], @(x) x(1));
assert(length(IVperm) <= N);

% Compute the number of points per vertex
SumOfWeightsPerVunorderd = accumarray(IV, W(:), [], @sum);
% Bring it in the right order
SumOfWeightsPerV = zeros(N, 1, class(SumOfWeightsPerVunorderd));
SumOfWeightsPerV(IVperm(IVperm ~= 0)) = SumOfWeightsPerVunorderd(IVperm ~= 0);

% Do the actual accumulation
for ii=1:size(X, 2)
  Yi = accumarray(IV, XW(:, ii), [], @sum);
  % Yi is still in the wrong oder and unreferenced vertices are not
  % included.
  % IVperm hold the right vertex indices.
  Y(IVperm(IVperm ~= 0), ii) = Yi(IVperm ~= 0);
end

% Normalize output
Y = Y ./ SumOfWeightsPerV;

% Set all outptu corresponding to unreferenced vertices to NaN
Y(SumOfWeightsPerV == 0, :) = NaN;

end

