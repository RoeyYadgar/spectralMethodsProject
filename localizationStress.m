function stressValue = localizationStress(patch)
patchDistanceMatrix = adjacency(patch,'weighted');
N = length(patchDistanceMatrix);
stressValue = 0;
for i = 1:N
    for j = (i+1):N
        stressValue = stressValue + (patchDistanceMatrix(i,j) - norm(patch.Nodes.Pos(i,:)-patch.Nodes.Pos(j,:))).^2;
    end
end
