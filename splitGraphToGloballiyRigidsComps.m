function patches = splitGraphToGloballiyRigidsComps(G)

N = size(G.Nodes,1);

edges = table2array(G.Edges);
patches = cell(N * 10,2);
patchNum = 1;
for i = 1:N
    nodes1 = edges(find(edges(:,1) == i),2);
    nodes2 = edges(find(edges(:,2) == i),1);
    neighboringNodes = union(nodes1,nodes2);
    nodeNames = cell(length(neighboringNodes)+1,1);
    nodeNames{1} = num2str(i);
    for j = 1:length(neighboringNodes)
        nodeNames{j+1} = num2str(neighboringNodes(j));
    end
    subStarGraph = subgraph(G,[i ;neighboringNodes]);
    subStarGraph.Nodes = nodeNames;
    
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