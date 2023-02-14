function [reflection,rotation,translation,residualNonReflected,residualReflected] = alignPatchesLSregis(patch1,patch2,commonNodes)
%computes rigid tranformation to align 2 patches that have 3 or more common
%nodes using least square registeration
%inputs:
%   patch1 - first patch
%   patch2 - second patch
%   commonNodes - array of common nodes between the two patches
%outputs:
%   reflection - reflection of the rigid tranformation (1 or -1)
%   rotation - rotation of the rigid transformation (e^(i*theta))
%   translation - translation of the rigid tranformation (complex number)

%we express each (x,y) points R2 as a complex number (x+iy), we can then
%reflect the graph (around the x-axis) by performing complex conjugate on
%all the data points and rotate the graph by multiplying by the points by
%e^(i*theta)

%Use the patches' subgraphs with their common nodes and transform their position to complex numbers:
[~,~,idx1] = intersect(commonNodes,patch1.Nodes.ID,'stable');
patch1Pos = patch1.Nodes.Pos(idx1);
[~,~,idx2] = intersect(commonNodes,patch2.Nodes.ID,'stable');
patch2Pos = patch2.Nodes.Pos(idx2);

%If the area of the convex hull of the patch is 0 we can't align the
%patches (in this case the function convhull will give an exception)
try
    [~,convhullArea] = convhull(real(patch1Pos),imag(patch1Pos));
    [~,convhullArea] = convhull(real(patch2Pos),imag(patch2Pos));
catch
    reflection = 0;
    rotation = 0;
    translation = 0;
    residualNonReflected = 0;
    residualReflected = 0;
    return 
end



%We solve 2 least squares problems, Ax = b, conj(A)x = b:
A = [patch2Pos ones(size(patch2Pos))];
b = patch1Pos;

x = ((A'*A)^-1)*(A'*b);
xResidual = norm(A*x-b);
A = conj(A);
y = ((A'*A)^-1)*(A'*b);
yResidual = norm(A*y-b);

if(xResidual < yResidual)  %Choose the relfection based on which least square solution has a lower residual error
    reflection = 1;
    result = x;
else
    reflection = -1;
    result = y;
end

%Normalize rotation to be of abs value of 1 (i.e. no rescaling) and fix the
%translation to account for the non rescaled position of patch2
rotation = result(1)/abs(result(1));
meanPatch2 = mean(patch2Pos);
if(reflection == -1)
    meanPatch2 = conj(meanPatch2);
end
translation = result(2) + (abs(result(1))-1) * meanPatch2 *rotation;


residualNonReflected = xResidual;
residualReflected = yResidual;



% plotPatch(patch1,patch1Pos)
% hold on
% plotPatch(patch2,patch2Pos*x(1)+x(2));
% plotPatch(patch2,conj(patch2Pos)*y(1)+y(2));
% axis('equal')
end
