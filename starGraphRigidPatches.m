function patches = starGraphRigidPatches(SG)
%Splits a star graph to maximally globally rigid components
%Input:
%   SG - star graph (assumed that first node is the star node)
%Outputs:
%   patches - cell array of patches (i.e. maximally globally rigid
%   components of the star graph)

%Remove all nodes of degree one
SG = subgraph(SG,find(degree(SG) > 1));

if(height(SG.Nodes) == 0) %If the graph is empty
    patches = {};
    return 
end

%Make a copy of the stargraph without the node and caluclate its
%bi-connected components
starNodeName = SG.Nodes{1,1};
SGcopy = subgraph(SG,[2:size(SG.Nodes,1)]);
bins2connected = biconncomp(SGcopy); %returns an index of the compontent for each edge in the graph
numComp = length(unique(bins2connected));

maxNodeSize = 32;
patches = cell(numComp*10,1);
counter = 1;
for i = 1:numComp
    patchesNodes = unique(SGcopy.Edges{bins2connected == i,1}(:)); %get all the nodes of the corresponding edges of the component 
    patch = subgraph(SG,[starNodeName ; patchesNodes]);
    if(height(patch.Nodes) >= maxNodeSize)
        nodeDeg = degree(patch);
        [~,maxDegNodes] = maxk(nodeDeg,maxNodeSize);
        patches{counter} = subgraph(patch,maxDegNodes);
        counter = counter + 1;
    else
        patches{counter} = patch;
        counter = counter+1;
    end
end

patches = patches(1:(counter-1));

end