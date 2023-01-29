function [patchReflection,patchRotation,patchTranslation] = generatePatchRelativeTransform(patches,graphDistance,rho)
%generates patch relative transformation matrix (for patches have
%enough nodes in common)
%inputs:
%   patches - N x 2 cell matrix of patches (first column) and their
%   corresponding star node (second column)
%   graphDistance - distance matrix of the original graph
%   rho - distance threshold of the Disc Graph Model
%outputs:
%   patchReflection - NxN matrix of relative reflective between patches
%   patchRotation - NxN matrix of relative rotation between patches
%   patchTranslation - NxN matrix of relative translation between patches



%Two patched can have common nodes only if their corresponding star nodes
%have atleast one common neighbor - ie if the matrix A^2 at the index i,j
%is greater than 0 (where A is the adjacency matrix of the original graph)
%This can speed up the process by elimanting the need to check for common nodes
%between most of the pairs of patches.
graphAdjacency = graphDistance > 0;
adjacencySquared = graphAdjacency^2;

N = size(patches,1); %number of total patches
patchReflection = sparse(N,N);
patchRotation = sparse(N,N);
patchTranslation = sparse(N,N);
C = sparse(N,N);
A = sparse(N,N);
Res = zeros(N,N,2);
Val = zeros(N,N,2);
for i = 1:N
    for j = (i+1):N
        if(adjacencySquared(patches{i,2},patches{j,2}) > 0)
            commonNodes = (intersect(patches{i,1}.Nodes.Name, patches{j,1}.Nodes.Name));
            C(i,j) = length(commonNodes);
            if(length(commonNodes) >= 3)
                [reflection,rotation,translation,Res(i,j,1),Res(i,j,2)] = alignPatchesLSregis(patches{i,1},patches{j,1},commonNodes);
                patchReflection(i,j) = reflection;
                patchRotation(i,j) = rotation;
                patchTranslation(i,j) = translation;
            elseif(length(commonNodes) == 2)
                [reflection,rotation,Val(i,j,1),Val(i,j,2)] = alignPatchesCombScore(patches{i,1},patches{j,1}, commonNodes,graphDistance,rho);
                patchReflection(i,j) = reflection;
                patchRotation(i,j) = rotation;
                
            end
            A(i,j) = length(commonNodes);
            A(j,i) = A(i,j);
        end
    end
end

patchReflection = patchReflection + patchReflection'; %Reflection between Pj and Pi is the same as Pi and Pj
patchRotation = patchRotation + patchRotation'; %Rotation between between Pj and Pi is conjugate of the rotation between Pi and Pj (the ' operator in matlab is conjugate transpose)
patchTranslation = patchTranslation - patchTranslation.'; %Translation between Pj and Pi is negative the tranlsation of Pi and Pj




end