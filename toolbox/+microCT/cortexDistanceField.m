function [cortexCenterDist, corticalThickness] = cortexDistanceField(perioMask, endoMask)
%CORTEXDISTANCEFIELD Computes the signed distance field from the cortex
%center.
%  [cortexCenterDist, corticalThickness] = cortexDistanceField(perioMask, endoMask) -
%    perioMask and endoMask muste be logical volumes representing the
%    periosteal and endosteal masks that are computed by `segmentCortex.
%    cortexCenterDist is the signed distance field representing the
%    distacne to the center of the cortex, negative values are inside,
%    positive outside.
%    corticalThickness represents the locl cortical thickness field. Values
%    outside the cortex are inf on the outside and -ing on the inside of
%    the bone.
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

pDist = bwdist(~perioMask);
eDist = bwdist(endoMask);

cortexCenterDist = eDist - pDist;
corticalThickness = pDist + eDist;
corticalThickness(~perioMask) = inf;
corticalThickness(endoMask) = -inf;

end

