function subpatches = splitPatch(patch)
%UNUSED
patchCopy = subgraph(patch,[2:size(patch.Nodes,1)]);

A = full(adjacency(patchCopy));
normalizedA = diag(sum(A,2).^-1)*A;
[EV,~] = eig(normalizedA);
secondEV = EV(:,2);
[partition,centers] = kmeans(secondEV,2);
threshold = mean(centers);

[~,boundaryNodes] = mink(abs(secondEV-threshold),2);


subpatches = {subgraph(patch,unique([1;find(partition == 1)+1;boundaryNodes+1])), ...
    subgraph(patch,unique([1;find(partition == 2)+1;boundaryNodes+1]))};

