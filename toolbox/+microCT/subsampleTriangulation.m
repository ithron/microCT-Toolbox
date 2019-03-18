function [Vnew, UVnew, Inew, Nnew] = subsampleTriangulation(V, F, N, varargin)
%SUBSAMPLETRIANGULATION Subsamples the given triangulation using uniform
%(non-random) barycentric subsampling.
%  [Vnew, UVnew, Inew] = subsampleTriangulation(V, F, N) - subsables the triangulation
%  given by vertices V and indices F. N is the number of subsamples per
%  edge.
%    Vnew is a Kx3 matrix containing the cartesian coordinates of the
%      samples
%    UVnew is a Kx2 matreix containing the barycnetric coordinates of the
%      samples
%    Inew is a K-vector containing the triangle indices of the samples
%
%  [Vnew, UVnew, Inew] = subsampleTriangulation(V, F, N, tolerance) - additionaly specify
%  a tolerance below that two points are considered equal, default to 1e-9;
%
%  [Vnew, UVnew, Inew, Nnew] = subsampleTriangulation(V, F, N, Normals, ...) -
%    also compute per-vertex normals and stores them in Nnew.
%
% The used algorithm is quite simple and not very efficient. Since it
% produces duplicate vertices at each edge, those duplicates are removed in
% a brute force filtering step afterwards. Not suited for large
% triangulations.
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

tol = 1e-9;
if not(isempty(varargin))
  if nargout == 3
    tol = varargin{1};
  else
    if numel(varargin) > 1
      tol = varargin{2};
    end
  end
end

if nargout == 4
  if isempty(varargin)
    error('Must specify per-vertex normals')
  else
    Normals = varargin{1};
    if size(Normals, 1) ~= size(V, 1)
      error('Normals must have the same size as V');
    end
  end
  
end

M = size(F, 1);

% Construct base matrix
VT = V(:, 1:3)';
T = [VT(:, F(:, 1)); VT(:, F(:, 2)); VT(:, F(:, 3))];
T = reshape(T, 3, []);

T = mat2cell(T, 3, 3 * ones(M, 1));
T = cell2mat(T');

if nargout == 4
  NT = Normals(:, 1:3)';
  TN = [NT(:, F(:, 1)); NT(:, F(:, 2)); NT(:, F(:, 3))];
  TN = reshape(TN, 3, []);
  
  TN = mat2cell(TN, 3, 3 * ones(M, 1));
  TN = cell2mat(TN');
end

% Construct barycentric coordinates
alpha = linspace(0, 1, N);
[a1, a2] = meshgrid(alpha, alpha);
a3 = 1 - a1 - a2;

valid = a3 >= -N^(-2);

a = [a1(valid)'; a2(valid)'; a3(valid)'];

% Compute new vertices
Vnew = T * a;
Vnew = reshape(Vnew, 3, []);

TI = reshape(repmat((1:M)', 1, size(a, 2)), [], 1);
UV = repmat(a(1:2, :)', M, 1);

if nargout == 4
  Nnew = TN * a;
  Nnew = reshape(Nnew, 3, []);
end

% Remove dubplicate vertices
[Vnew, IA] = microCT.filterClosePointsQuick(Vnew', tol);

Inew = TI(IA);
UVnew = UV(IA, :);

if nargout == 4
  
  Nnew = Nnew(:, IA)';
  
end

end

