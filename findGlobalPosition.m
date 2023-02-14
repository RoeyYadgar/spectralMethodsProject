function posNodes = findGlobalPosition(M,patches,reflections,rotations)
%finds the global position of the nodes (up to a rigid transformation)
%inputs:
%   M - number of total Nodes in the data
%   patchReflection - NxN matrix of relative reflection between patches
%   patchRotation - NxN matrix of relative rotation between patches
%   patchTranslation - NxN matrix of relative translation between patches
%outputs:
%   posNodes - Mx2 matrix of x,y coordinates of the M nodes


edgeMat = sparse(M,M);
edgeCounter = sparse(M,M);
N = size(reflections,1);

for i = 1:N
    if(reflections(i) ~= 0) %patches that aren't in the biggest connected component (of the patch graph) cannot be aligned and will have reflection,rotation = 0 
        patchEdges = adjacency(patches{i,1});
        patchPos = patches{i,1}.Nodes.Pos;
        patchID = patches{i,1}.Nodes.ID;
        %Align the patch to the global coordinate system by applying it's
        %inverse transformation
        if(reflections(i) == -1)
            patchPos = conj(patchPos);
        end
        patchPos = patchPos * conj(rotations(i));

        %Loop over all the edges of the patch and sum the edge 'vector'
        %(complex valued) 
        for j = 1:(size(patchEdges,1))
            for k = (j+1):size(patchEdges,2)
                if(patchEdges(j,k) == 1)
                    node1 = patchID(j);
                    node2 = patchID(k);
                    edgeMat(node1,node2) = edgeMat(node1,node2) + patchPos(j)-patchPos(k);
                    edgeCounter(node1,node2) = edgeCounter(node1,node2) +1;
                end
            end
        end
    end
end

%while summing the same edge over different patches we might sum them in the
%opposite direction (ie patchPos(j) - patchPos(i) instead of patchPos(i) -
%patchPos(j)) and we will store it in edgeMat(j,i) instead of edgeMat(i,j)
% - to get the total sum with the same direction we can substract the
% transposed matrix from itself. similarly edge counter can be summed with
% its transposed matrix.
%we can then take the upper triangular part of it since the j-th,i-th value will
%generate the same equation as the i-th,j-th value in the Least Squares
%problem (with a different sign)
edgeMat = triu(edgeMat - edgeMat.');
edgeCounter = triu(edgeCounter + edgeCounter');

%Build the Least Squares Matrix T - each row corresponds to an edge in the
%graph 
[a,b] = find(edgeCounter > 0);
edgeVector = zeros(size(a));
T = sparse(size(a,1),M);
for i = 1:length(a)
    edgeVector(i) = edgeMat(a(i),b(i));
    T(i,a(i)) = edgeCounter(a(i),b(i));
    T(i,b(i)) = -edgeCounter(a(i),b(i));
end

%Solve the Least Squares problem 
posNodesComplex = T\edgeVector;
posNodes = zeros(M,2);
posNodes(:,1) = real(posNodesComplex);
posNodes(:,2) = imag(posNodesComplex);



end