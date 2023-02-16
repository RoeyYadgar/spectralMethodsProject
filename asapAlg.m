function [aneScore, alignmentStats] = asapAlg(dataPoints,rho,eta)

% profile on
G = generateGraphDiscModel(dataPoints,rho);
W = sparse(adjacency(G,'weighted'));
patches = splitGraphToGloballiyRigidsComps(G);

[patchReflection,patchRotation,~] = generatePatchRelativeTransform(patches,W,rho);
[reflections,rotations,correctedPatchRotation] = findGlobalTransformation(patchReflection,patchRotation);

[posNodes,nonLocalizedNodes] = findGlobalPosition(length(dataPoints),patches,reflections,rotations);
aneScore = ANE(dataPoints,posNodes,nonLocalizedNodes);

[refsGT, rotsGT] = transformationGT(patches,dataPoints);
alignmentStats = alignmentStatistics(reflections,rotations,patchReflection,correctedPatchRotation,refsGT,rotsGT);

% profile viewer


end