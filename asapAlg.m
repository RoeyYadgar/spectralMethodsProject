function [aneScore, alignmentStats] = asapAlg(dataPoints,rho,eta)
%Solves the SNL problem using the ASAP algorithm
%inputs:
%   dataPoints - Nx2 matrix with position of the data 
%   rho - distance threshold constant
%   eta - multiplictive noise level
%outputs:
%   aneScore - ANE of the estimated position of the nodes
%   alignmentStats - 7x1 cell array of alignment statistics from table 4, page
%       33 in paper

%Generate the graph
G = generateGraphDiscModel(dataPoints,rho,eta);
W = sparse(adjacency(G,'weighted'));

%Split the graph into globally rigid patches and localize each patch
patches = splitGraphToGloballiyRigidsComps(G);

%Align each pair of patches with common nodes and find a global alignment of reflections and rotations
[patchReflection,patchRotation,~] = generatePatchRelativeTransform(patches,W,rho);
[reflections,rotations,correctedPatchRotation] = findGlobalTransformation(patchReflection,patchRotation);

%Find global position of the nodes (up to some rigid transformation)
[posNodes,nonLocalizedNodes] = findGlobalPosition(length(dataPoints),patches,reflections,rotations);

%Compute aneScore and alignmentStats
aneScore = ANE(dataPoints,posNodes,nonLocalizedNodes);
[refsGT, rotsGT] = transformationGT(patches,dataPoints);
alignmentStats = alignmentStatistics(reflections,rotations,patchReflection,correctedPatchRotation,refsGT,rotsGT);


end
