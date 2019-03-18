function E = meshEdges(F)
%MESHEDGES Output unique edges of a mesh given py its index set
%
%  E = meshEdges(F) - returns unique edges of mesh defined by the face
%  matrix F.
%    F is a Nx3 matrix of indices
%    E is a Mx2 matrix of indices
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

E = reshape([F, F(:, 1)]', 2, []);
E = sort(E, 1);
E = unique(E', 'rows');

end

