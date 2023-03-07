function patches = splitGraphToGloballiyRigidsComps(G)
%splits the orignal graph into globally rigid patches and localizes each
%patch
%inputs:
%   G - original graph (with the assumption of the disc graph model)
%outputs:
%   patches - Nx2 cell matrix of patches (first column) and their
%       corresponding star node (second column)

N = size(G.Nodes,1);
W = adjacency(G,'weighted');
patches = cell(N * 10,2);
patchNum = 1;
for i = 1:N
    %get 1-hop neighborhood star graph
    neighboringNodes = find(W(i,:) > 0)';
    nodeNames = cell(length(neighboringNodes)+1,1);
    nodeNames{1} = num2str(i);
    for j = 1:length(neighboringNodes)
        nodeNames{j+1} = num2str(neighboringNodes(j));
    end
    subStarGraph = subgraph(G,[i ;neighboringNodes]);
    subStarGraph.Nodes = nodeNames;
    
    %split patch to globally rigid components and localize each patch
    starGraphPatches = starGraphRigidPatches(subStarGraph);
    patches(patchNum : (patchNum + length(starGraphPatches)-1),1) = starGraphPatches(:);
    patches(patchNum : (patchNum + length(starGraphPatches)-1),2) = {i};
    for p = 1:length(starGraphPatches)
        patches{patchNum + p - 1,1}.Nodes.Pos = localizePatch(starGraphPatches{p});
        patches{patchNum + p - 1,1}.Nodes.ID = cellfun(@str2num,starGraphPatches{p}.Nodes.Name);     
    end
    patchNum = patchNum + length(starGraphPatches);
    
    
end
patches = patches(1:(patchNum-1),:);




end