function estimationError = ANE(realPos,estimatedPos,nonLocalizedNodes)
%

localizedNodes = true(length(realPos),1);
localizedNodes(nonLocalizedNodes) = false;
realPos = realPos(localizedNodes,:);
estimatedPos = estimatedPos(localizedNodes,:);

N = length(realPos);
graph1 = graph(zeros(N,N));
graph1.Nodes.Pos = realPos(:,1) + 1i*realPos(:,2);
graph1.Nodes.ID = (1:N)';

graph2 = graph1;
estimatedPosComplex = estimatedPos(:,1) + 1i*estimatedPos(:,2);
graph2.Nodes.Pos = estimatedPosComplex;

[reflection,rotation,translation] = alignPatchesLSregis(graph1,graph2,(1:N)');

if(reflection == -1)
    estimatedPosComplex = conj(estimatedPosComplex);
end
estimatedPosComplex = rotation*estimatedPosComplex + translation;

estimatedPos = [real(estimatedPosComplex) imag(estimatedPosComplex)];

meanPos = mean(realPos);


estimationError = sqrt(sum(vecnorm(realPos-estimatedPos,2,2).^2))/sqrt(sum(vecnorm(realPos-meanPos,2,2).^2));

plotPatch(graph1)
hold on
plotPatch(graph2,estimatedPos)
end
