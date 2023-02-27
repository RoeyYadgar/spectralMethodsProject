function nodesPos = localizePatch(patch)
%calculates position of nodes of the patch 
%inputs:
%   patch - a patch graph (small globally rigid graph computed in the
%   function starGraphRigidPatches)
%ouputs:    
%   nodesPos - N x 2 matrix containing the position of the nodes in the
%   graph 

N = height(patch.Nodes);
distMatrix = full(adjacency(patch,'weighted')); %The weight adjacency matrix contains the distances between the nodes
estimatedDistMatrix = distMatrix;

%Estimate distance between nodes where distance measurement is missing
for i = 1:N
    for j  = (i+1):N
        if(distMatrix(i,j) == 0)
            ikDist = distMatrix(i,:); %distance between the node i to the rest of the nodes
            jkDist = distMatrix(j,:); %distance between the node j o the rest of the nodes
            
            %Calculate Lower and Upper bound on the missing distance
            %measurement:
            %lower bound is the maximal distance beween the nodes i and j
            %to the rest of the nodes (bound is only vaid because of
            %the disc graph model assumption)
            distLowerBound = max([ikDist jkDist]); 
            %upper bound is the minimal d_ik + d_kj where k is a node that
            %is connected to both i and j (triangle inequality)
            ikDist(ikDist==0) = inf; %if there is no edge between i and k set the distance to inf
            jkDist(jkDist==0) = inf; %if there is no edge between j and k set the distance to inf
            distUpperBound = min(ikDist+jkDist);
            estimatedDistMatrix(i,j) = 0.5*(distLowerBound + distUpperBound);
            estimatedDistMatrix(j,i) = estimatedDistMatrix(i,j);
        end
        
    end
end

%Perfrom classical MDS:
nodesPos = cmdscale(estimatedDistMatrix,2);

if(size(nodesPos,2) ~= 2) %sometime a graph embedding can be on a single line (in this case dont perform majorization techinque)
    nodesPos = nodesPos(:,1) + 1i * 0;
    return
end

%Perfrom iteartive majorization techinque to minimize the stress function:
nodeDegs = degree(patch);
updatedNodesPos = zeros(size(nodesPos));
iterativeDifference = 1;
while(iterativeDifference > 1e-5) %perfrom iterations until the chance between consectuve iterations is small 
    for i = 1:N
        newNodePosSum = [0 0];
        for j = 1:N
           if(distMatrix(i,j) > 0)
               Pi = nodesPos(i,:);
               Pj = nodesPos(j,:);
               Dij = distMatrix(i,j);
               normPiPj = norm(Pi-Pj);
               if(normPiPj ~= 0) %inverse should be 0 if the norm itself is 0
                   invNormPiPj = 1/normPiPj;
               else
                   invNormPiPj = 0;
               end 
               newNodePosSum = newNodePosSum + (Pj + Dij*(Pi-Pj)*invNormPiPj);  
           end
        end
        updatedNodesPos(i,:) = newNodePosSum / nodeDegs(i);
        
    end
    
    iterativeDifference = mean(vecnorm((nodesPos-updatedNodesPos)'));
    nodesPos = updatedNodesPos;
    
end

%Convert nodePos to complex pos
nodesPos = nodesPos(:,1) + 1i * nodesPos(:,2);


end