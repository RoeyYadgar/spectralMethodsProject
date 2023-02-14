function [reflections,rotations] = findGlobalTransformation(patchReflection,patchRotation,Adj)

N = length(patchReflection);
z = patchReflection;
connectedPatches = find(conncomp(graph(abs(z))) == 1);
z = z(connectedPatches,connectedPatches);
Z = sparse(diag(sum(abs(z),2).^(-1)))*z;
[Vref,~] = eigs(Z,1);
reflections = zeros(N,1);
reflections(connectedPatches) = sign(Vref);
plot(abs(Vref));

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


x = 1;
