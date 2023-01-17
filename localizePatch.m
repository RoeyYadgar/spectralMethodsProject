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
nodesPos = mdscale(estimatedDistMatrix,2);



end