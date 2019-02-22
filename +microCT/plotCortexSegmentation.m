function f = plotCortexSegmentation(imageSlice, perioMaskSlice, endoMaskSlice)
% PLOTCORTEXSEGMENTATION Visualizes a cortex segmentation
%  plotCortexSegmentation(imageSlice, perioMaskSlice, endoMaskSLice) -
%    overlays the image slice with the periosteal mask perioMaskSlice and the
%    endosteal mask endoMaskSlice.
%
%  f = plotCortexSegmentation(imageSlice, perioMaskSlice, endoMaskSLice) -
%    Returns the figure handle f
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.

imageSlice = single(imageSlice - min(imageSlice(:))) / ...
  single(max(imageSlice(:)) - min(imageSlice(:)));

rFactor = single(perioMaskSlice);
gFactor = single(endoMaskSlice);
bFactor = single(~perioMaskSlice & ~endoMaskSlice);

colImg = ...
  cat(3, imageSlice .* rFactor, imageSlice .* gFactor, imageSlice .* bFactor);

if nargout > 0
  f = imagesc(colImg);
else
  imagesc(colImg);
end

axis equal

end
