function estimationError = ANE(realPos,estimatedPos,nonLocalizedNodes)
%computes ANE score of the estimated position of the nodes
%inputs:
%   realPos - Nx2 ground truth position of the nodes
%   estimatedPos - Nx2 estimated position of the nodes
%   nonLocalizedNodes - array of indices of nodes that were not localized
%outputs:
%   estimationError - ANE score of the estimated position

%Take only the nodes that were localized
localizedNodes = true(length(realPos),1);
localizedNodes(nonLocalizedNodes) = false;
realPos = realPos(localizedNodes,:);
estimatedPos = estimatedPos(localizedNodes,:);

%Generate graphs with and set their position to the ground truth and the
%estimaed position
N = length(realPos);
graph1 = graph(zeros(N,N));
graph1.Nodes.Pos = realPos(:,1) + 1i*realPos(:,2);
graph1.Nodes.ID = (1:N)';

graph2 = graph1;
estimatedPosComplex = estimatedPos(:,1) + 1i*estimatedPos(:,2);
graph2.Nodes.Pos = estimatedPosComplex;

%Use least squares method to align the ground truth and the estimated
%position
[reflection,rotation,translation] = alignPatchesLSregis(graph1,graph2,(1:N)');

%apply the transformation to align
if(reflection == -1)
    estimatedPosComplex = conj(estimatedPosComplex);
end
estimatedPosComplex = rotation*estimatedPosComplex + translation;

%Compute ANE
estimatedPos = [real(estimatedPosComplex) imag(estimatedPosComplex)];
meanPos = mean(realPos);
estimationError = sqrt(sum(vecnorm(realPos-estimatedPos,2,2).^2))/sqrt(sum(vecnorm(realPos-meanPos,2,2).^2));

plotPatch(graph1)
hold on
plotPatch(graph2,estimatedPos)
end
