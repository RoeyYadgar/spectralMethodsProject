function [reflections,rotations] = findGlobalTransformation(patchReflection,patchRotation)
%finds the global transformation of the patches (up to a rigid transformation)
%inputs:
%   patchReflection - NxN matrix of relative reflection between patches
%   patchRotation - NxN matrix of relative rotation between patches
%outputs:
%   reflections - Nx1 vector of global reflection of each patch
%   rotations - Nx1 vector of global rotation of each patch


N = length(patchReflection);
z = patchReflection;
connectedPatches = find(conncomp(graph(abs(z))) == 1);
z = z(connectedPatches,connectedPatches);
Z = sparse(diag(sum(abs(z),2).^(-1)))*z;
[Vref,~] = eigs(Z,1);
reflections = zeros(N,1);
reflections(connectedPatches) = sign(Vref);

r = triu(patchRotation);
for i = 1:N
    if(reflections(i) == -1)
        r(i,:) = conj(r(i,:));
    end
end
r = r(connectedPatches,connectedPatches);
r = r+r';
R = sparse(diag(sum(abs(r),2).^(-1)))*r;
[Vrot,~] = eigs(R,1);
rotations = zeros(size(reflections));
rotations(connectedPatches) = Vrot(:,1)./abs(Vrot(:,1));


