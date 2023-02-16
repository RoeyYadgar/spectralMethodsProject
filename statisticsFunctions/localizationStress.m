function stressValue = localizationStress(patch)
patchDistanceMatrix = adjacency(patch,'weighted');
N = length(patchDistanceMatrix);
E = 0;
stressValue = 0;
for i = 1:N
    for j = (i+1):N
        if(patchDistanceMatrix(i,j) > 0)
            stressValue = stressValue + (patchDistanceMatrix(i,j) - norm(patch.Nodes.Pos(i,:)-patch.Nodes.Pos(j,:))).^2;
            E = E+1;
        end
    end
end
stressValue = stressValue/E;
