function splitGraphToGloballiyRigidsComps(G)

N = size(G.Nodes,1);

edges = table2array(G.Edges);
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
    
    
end




end