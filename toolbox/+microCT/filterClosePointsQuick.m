function [Pf, IA] = filterClosePointsQuick(P, r)
%FILTERCLOSEPOINTS Filter out points that have a distance less that r
%inbetween.
%  [Pf, IA] = filterClosePoints(P, r)
%    r is the filter radius
%    P is a NxM matrix of M-dimensional points
%    Pf is a KxM matrix of the remaining K M-dimensional points.
%    IA is a index matrix, so that Pf = P(IA, :)
%
% It is not guaranteed that this function filters out all points that meet
% the filter condition, but it is fast. For exact filtering use
% uniquetol().
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

  function I = medPointIdx(X)
    PX = P(X, :);
    centroid = mean(PX, 1);
    dist = sum((PX - centroid).^2, 2);
    [~, iMin] = min(dist, [], 1);
    I = X(iMin);
  end

N = size(P, 1);
griddedP = round(P ./ (sqrt(3) * r));

[~, ~, IB] = unique(griddedP, 'rows');

% For each gridded point, find the central point of all points inside the
% grid cell
IA = accumarray(IB, (1:N)', [], @medPointIdx);
Pf = P(IA, :);

end

