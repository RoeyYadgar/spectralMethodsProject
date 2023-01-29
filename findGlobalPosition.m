function posNodes = findGlobalPosition(M,patches,reflections,rotations)
%finds the global position of the nodes (up to a rigid transformation)
%inputs:
%   M - number of total Nodes in the data
%   patchReflection - NxN matrix of relative reflection between patches
%   patchRotation - NxN matrix of relative rotation between patches
%   patchTranslation - NxN matrix of relative translation between patches
%outputs:
%   posNodes - Mx2 matrix of x,y coordinates of the M nodes
% % %   reflections - Nx1 vector of global reflection of each patch
% % %   rotations - Nx1 vector of global rotations of each patch
% % %   translations - Nx1 vector of global translations of each patch

% %for relfection and rotation matrices - normalize each row by the degree of
% %the patch and find biggest eigenvector
% normalizationMat = diag(sum(abs(patchReflection),2).^-1);
% [evRef,~] = eigs(normalizationMat*patchReflection,1);
% reflections = sign(evRef);
% [evRot,~] = eigs(normalizationMat*patchRotation,1);
% rotations = evRot./abs(evRot);

edgeMat = sparse(M,M);
edgeCounter = sparse(M,M);
N = size(reflections,1);

for i = 1:N
    patchEdges = adjacency(patches{i,1});
    patchPos = patches{i,1}.Nodes.Pos;
    patchPos = patchPos(:,1) + 1i*patchPos(:,2);
    if(reflections(i) == -1)
        patchPos = conj(patchPos);
    end
    patchPos = patchPos * rotations(i);
    for j = 1:(size(patchEdges,1))
        for k = (j+1):size(patchEdges,2)
            if(patchEdges(j,k) == 1)
                node1 = str2double(patches{i,1}.Nodes.Name{j});
                node2 = str2double(patches{i,1}.Nodes.Name{k});
                edgeMat(node1,node2) = edgeMat(node1,node2) + patchPos(j)-patchPos(k);
                edgeCounter(node1,node2) = edgeCounter(node1,node2) +1;
            end
            
            
        end
        
        
    
    
    
    end
end

edgeMat = triu(edgeMat - edgeMat.');
edgeCounter = triu(edgeCounter + edgeCounter');

[a,b] = find(edgeCounter > 0);
edgeVector = zeros(size(a));
T = sparse(size(a,1),M);
for i = 1:length(a)
    edgeVector(i) = edgeMat(a(i),b(i));
    T(i,a(i)) = edgeCounter(a(i),b(i));
    T(i,b(i)) = -edgeCounter(a(i),b(i));
end

posNodesComplex = T\edgeVector;
posNodes = zeros(M,2);
posNodes(:,1) = real(posNodesComplex);
posNodes(:,2) = imag(posNodesComplex);



end