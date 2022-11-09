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
nBytes = (isa(cortexCenterDist, 'single')) * 4 + (isa(cortexCenterDist, 'double')) * 8;
% Scale voluem if too large (there seems to be a 2GB limit)
maxN = 2*1024^3 / nBytes;
scale = 1;
if nVox > maxN
  scale = nthroot(maxN / nVox, 3);
  scale = min(floor(size(cortexCenterDist) * scale) ./ size(cortexCenterDist));
  
  cortexCenterDist = imresize3(cortexCenterDist, scale, 'linear');
end


[rootDir, ~, ~] = fileparts(mfilename('fullpath'));
oldPath = addpath(sprintf('%s/iso2mesh/', rootDir));

% Mesh volume
r = max(size(cortexCenterDist))^3 / 2;
c = size(cortexCenterDist) / 2;
[V, F] = vol2restrictedtri(cortexCenterDist, -1e-9, c, r, 30, 5, 0.1, 500000);
V = single(V);
if length(V) <= intmax('uint16')
  F = uint16(F);
elseif length(V) <= intmax('uint32')
  F = uint32(F);
else
  F = uint64(F);
end


% Undo scaling
V = V / scale + 0.5;
V = [V(:, 2), V(:, 1), V(:, 3)];

[V, F] = meshcheckrepair(V, F, 'dup');
[V, F] = meshcheckrepair(V, F, 'isolated');
[V, F] = meshcheckrepair(V, F, 'deep');
[V, F] = meshcheckrepair(V, F, 'meshfix');

path(oldPath);

% Sample corticalThickness
C = interp3(corticalThickness, V(:, 1), V(:, 2), V(:, 3), 'linear');

end
