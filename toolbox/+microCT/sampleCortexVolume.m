function [CtBMD, Ns] = sampleCortexVolume(Img, V, N, CtTh)
%SAMPLECORTEXVOLUME Samples the cortex volume along profiles through the
%given vertices along the given normals and return the BMD values at the
%sampling positions
%
%  [CtBMD, Ns] = sampleCortexVolume(Img, V, N, CtTh)
%    Img - volume to sample (grayscale floating point)
%    V - Nx3 matrix containing vertex positions
%    N - Nx3 matrix of per-vertex normals
%    CtTh - N21 vector of per-vertex cortical thickness values
%
%    CtBMD - Nx1 vector of mean BMD values accumulated along profiles through V
%      along N.
%    Ns - Nx1 vector containing the number of samples acquired per profile
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

CtTh = ceil(CtTh/2);

maxTh = max(CtTh(isfinite(CtTh)));
minTh = -maxTh;

dTh = 0.5;

ths = [minTh : dTh : maxTh, 0];

V = single(V);
N = single(N);
CtTh = single(CtTh);
CtTh(~isfinite(CtTh)) = NaN;

if length(ths) < intmax('uint8')
  Nclass = 'uint8';
elseif length(ths) < intmax('uint16')
  Nclass = 'uint16';
elseif length(ths) < intmax('uint32')
  Nclass = 'uint32';
else
  Nclass = 'uint64';
end

CtBMD = zeros(size(V, 1), 1, 'double');
Ns = zeros(size(V, 1), 1, Nclass);

for ii=1:length(ths)
  
  thi = ths(ii);
  
  valid = abs(CtTh) >= abs(thi);
  
  Vi = V(valid, :) + thi * N(valid, :);
  
  Ni = ones(size(Vi, 1), 1, Nclass);
  CtBMDi = interp3(single(Img), Vi(:, 1), Vi(:, 2), Vi(:, 3), 'nearest');
  valid2 = isfinite(CtBMDi);
  CtBMDi(~valid2) = 0;
  Ni(~valid2) = 0;
  
  Ns(valid) = Ns(valid) + Ni;
  
  CtBMD(valid) = CtBMD(valid) + CtBMDi;
  
  fprintf('%.1f%%\n', 100 * ii / length(ths));
  
end

CtBMD = single(CtBMD ./ double(Ns));

end

