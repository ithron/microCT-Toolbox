function [F, V, C] = meshCortex(cortexCenterDist, corticalThickness, varargin)
%MESHCORTEX Meshes the surface of the cortex center, colored with the lcoal
%cortical thickness
%  [F, V, C] = meshCortex(cortexCenterDist, corticalThickness) - Outputs
%    the triangle index list F (Mx3), the vertex list V (Nx3) and the
%    per-vertex cortical thickness (color) values C (Nx1).
%    cortexCenterDist and corticalThickness are the output of
%    cortexDistanceField.
%
%  [F, V, C] = meshCortex(cortexCenterDist, corticalThickness, smoothing) -
%    smooth the distance field before the meshing to avoid aliasing.
%    smoothing is the sigma of a gaussian filter kernel.
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.


if nargin == 3
  s = varargin{1};
  kS = floor((3 * s) / 2) * 2 + 1;
  cortexCenterDist = smooth3(cortexCenterDist, 'gaussian', [1 1 1] * kS, s);
end


nVox = numel(cortexCenterDist);

% Scale voluem if too large
scale = 1;
if nVox > 1e9
  scale = 1e9 / nVox;
  
  cortexCenterDist = imresize3(cortexCenterDist, scale);
end


[rootDir, ~, ~] = fileparts(mfilename('fullpath'));
oldPath = addpath(sprintf('%s/iso2mesh/', rootDir));

% Mesh volume
r = max(size(cortexCenterDist))^3 / 2;
c = size(cortexCenterDist) / 2;
[V, F] = vol2restrictedtri(double(cortexCenterDist), -1e-9, c, r, 30, 2, 2, 50000);

path(oldPath);

% Undo scaling
V = (V + 1) / scale - 0.5;
V = [V(:, 2), V(:, 1), V(:, 3)];

% Sample corticalThickness
C = interp3(corticalThickness, V(:, 1), V(:, 2), V(:, 3), 'cubic');

end