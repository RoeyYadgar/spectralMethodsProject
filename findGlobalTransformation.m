function [reflections,rotations,correctedPatchRotation] = findGlobalTransformation(patchReflection,patchRotation)
%finds the global transformation of the patches (up to a rigid transformation)
%inputs:
%   patchReflection - NxN matrix of relative reflection between patches
%   patchRotation - NxN matrix of relative rotation between patches
%outputs:
%   reflections - Nx1 vector of global reflection of each patch
%   rotations - Nx1 vector of global rotation of each patch
%   correctedPatchRotation - NxN patchRotation matrix after reflection
%       correction

%the patch graph is not necessarily connected (usually there will be one
%large connected component and a few very small connected component) - we
%therefore take the biggest connected component and find the global
%tranformation of each patch in that component
N = length(patchReflection);
z = patchReflection;
connectedPatches = find(conncomp(graph(abs(z))) == 1);
z = z(connectedPatches,connectedPatches);

Z = sparse(diag(sum(abs(z),2).^(-1)))*z;
[Vref,~] = eigs(Z,1);
reflections = zeros(N,1);
reflections(connectedPatches) = sign(Vref);

%while measuring the relative transformation between patches, the rotation
%between the p1 and p2 is the inverse rotation between p1_reflected and p2.
%to fix it we update the matrix patchRotation by applying conjugate to
%every row where the global reflection of the corresponding patch is -1
r = triu(patchRotation);
for i = 1:N
    if(reflections(i) == -1)
        r(i,:) = conj(r(i,:));
    end
end
r = r(connectedPatches,connectedPatches);
r = r+r';
correctedPatchRotation = zeros(length(reflections),length(reflections));
correctedPatchRotation(connectedPatches,connectedPatches) = r;

R = sparse(diag(sum(abs(r),2).^(-1)))*r;
[Vrot,~] = eigs(R,1);
rotations = zeros(size(reflections));
rotations(connectedPatches) = Vrot(:,1)./abs(Vrot(:,1));

% %Generate plot of eigenvectors of reflections and rotations
% figure
% Vref = Vref/sum(abs(Vref))*length(Vref);
% histogram(Vref,10);
% title('histogram of entries of the top eigenvector of the matrix Z')
% figure
% Vrot = Vrot/sum(abs(Vrot))*length(Vrot);
% scatter(real(Vrot),imag(Vrot));
% title('scatter plot of entries of the top eigenvector of the matrix R')


