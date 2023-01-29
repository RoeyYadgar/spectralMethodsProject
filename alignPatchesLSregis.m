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
patch1 = subgraph(patch1,commonNodes);
patch1Pos = patch1.Nodes.Pos(:,1) + 1i *  patch1.Nodes.Pos(:,2);
patch2 = subgraph(patch2,commonNodes);
patch2Pos = patch2.Nodes.Pos(:,1) + 1i *  patch2.Nodes.Pos(:,2);


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
end
