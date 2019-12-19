function varargout = segmentCortex(img, th1, th2, smoothing, varargin)
% SEGMENTCORTEX Computes a cortex segmentation mask
%    [perioMask, endoMask] = segmentCortex(img, th1, th2, smoothing) -
%      Computes a periosteal mask (perioMask) that includes everything inside
%      the periosteal surface and a endosteal mask (endoMask) that includes
%      everything inside the endosteal surface. The cortex mask can be
%      obtaine with `periomask & ~endoMask`.
%      th1 is the threshold used for the periosteal mask, th2 is the
%      threshold for the endosteal mask.
%      smoothing is a smoothign factor applied to the endosteal mask.
%
%    cortexMask = segmentCortex(img, th1, th2, smoothing) - only computes
%      the cortex mask
%
%    [perioMask, endoMask] = segmentCortex(img, th1, th2, smoothing, type)
%    cortexMask = segmentCortex(img, th1, th2, smoothing, type)
%      -- same as above, but scanenr tpye can be selected. `microCT` (default)
%      micro CT configuration, `XCT` XtremeCT (or HR-pQCT) scanner
%      configuation.
%
% This file is part of the 'microCT-Toolbox' project.
% Author: Stefan Reinhold
% Copyright: Copyright (C) 2019 Stefan Reinhold  -- All Rights Reserved.
%            You may use, distribute and modify this code under the terms of
%            the AFL 3.0 license; see LICENSE for full license details.


if not(isempty(varargin))
  if not(strcmpi(varargin{1}, 'microCT') || strcmpi(varargin{1}, 'XCT'))
    error(['Invalid type: "' varargin{1} '"']);
  else
    type = varargin{1};
  end
else
  type = 'microCT';
end

if strcmpi(type, 'XCT')
  padding = 13;
  medFiltSize = 3;
  se1Radius = 11;
  se2Radius = 5;
elseif strcmpi(type, 'microCT')
  padding = 32;
  medFiltSize = 7;
  se1Radius = 30;
  se2Radius = 10;
end

paddedImg = padarray(img, padding * ones(1, 3));

se1 = strel('sphere', se1Radius);
se2 = strel('sphere', se2Radius);

% Compute periosteal mask
perioMask = paddedImg > th1;
perioMask = medfilt3(perioMask, medFiltSize * ones(1, 3));
perioMask = imdilate(perioMask, se1);
CC1 = bwconncomp(~perioMask);
perioMask = false(size(perioMask));
perioMask(CC1.PixelIdxList{1}) = true;
perioMask = ~perioMask;
perioMask = imerode(perioMask, se1);

% Compute endosteal mask
endoMask = paddedImg > th2;
endoMask = ~perioMask | endoMask;
endoMask = imgaussfilt3(single(endoMask), smoothing) > 0.5;
endoMask = imdilate(~endoMask, se2);
CC2 = bwconncomp(~endoMask);
endoMask = false(size(endoMask));
endoMask(CC2.PixelIdxList{1}) = true;
endoMask = ~endoMask;
endoMask = imerode(endoMask, se2);

% undo padding

perioMask = perioMask(...
  padding + 1 : end - padding, ...
  padding + 1 : end - padding, ...
  padding + 1 : end - padding);

endoMask = endoMask(...
  padding + 1 : end - padding, ...
  padding + 1 : end - padding, ...
  padding + 1 : end - padding);

if nargout == 2
  varargout{1} = perioMask;
  varargout{2} = endoMask;
elseif nargout == 1
  varargout{1} = perioMask & ~endoMask;
else
  error('Unsupported number of outputs');
end

end

