function [CtTh, CtBMD] = perVertexCtThCtBMD(F, V, C, Img, CtThField)
%PERVERTEXCTTHCTBMD Summary of this function goes here
%   Detailed explanation goes here

[Vnew, UVnew, Inew, Nnew, Cnew] = microCT.genMeshSamples(F, V, C);

Img = single(Img);
Img(~isfinite(CtThField)) = NaN;

[CtBMD, ~]  = microCT.sampleCortexVolume(Img, Vnew, Nnew, Cnew);

Y = microCT.accumVertexData(F, UVnew, Inew, [CtBMD, Cnew]);

CtTh = Y(:, 2);
CtBMD = Y(:, 1);

end

