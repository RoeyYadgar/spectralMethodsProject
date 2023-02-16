function statistics = alignmentStatistics(reflections,rotations,patchReflection,patchRotation,reflectionsGT,rotationsGT)
%Calculates alignment statistics (table 4 in page 33 of the paper)

refErr = reflections.*reflectionsGT;
tau = sum(refErr == -1)/sum(refErr ~= 0);

%reflections vector idealy should be either reflectionsGT or -reflectionsGT (in the case of -reflectionsGT tau will be 1-errorRate)
if(tau > 1 - tau)
    tau = 1 - tau;
    reflections = -reflections;
    rotations = conj(rotations);
    patchRotation = conj(patchRotation);
end
tau = 100*tau;


patchReflectionGT = sparse(diag(reflectionsGT)*abs(patchReflection)*diag(reflectionsGT));
patchRefErr = patchReflectionGT.*patchReflection;
epsZinp = 100*sum(sum(patchRefErr == -1))/sum(sum(patchRefErr ~= 0));

patchReflectionEig = sparse(diag(reflections)*abs(patchReflection)*diag(reflections));
patchRefErrEig = patchReflectionGT.*patchReflectionEig;
epsZeig = 100*sum(sum(patchRefErrEig == -1))/sum(sum(patchRefErrEig ~= 0));

patchRotationGT = sparse(diag(rotationsGT)*abs(patchRotation)*diag(conj(rotationsGT)));
patchRotErr = patchRotationGT.*conj(patchRotation);
epsRinp = mean(abs(angle(patchRotErr(patchRotErr ~= 0))))/pi*180;
outliersRinp = 100*sum(abs(angle(patchRotErr(patchRotErr ~= 0))/pi*180) > 10) / sum(sum(patchRotErr ~= 0));

patchRotationEig = sparse(diag(rotations)*abs(patchRotation)*diag(conj(rotations)));
patchRotErrEig = patchRotationGT.*conj(patchRotationEig);
epsReig = mean(abs(angle(patchRotErrEig(patchRotErrEig ~= 0))))/pi*180;
outliersReig = 100*sum(abs(angle(patchRotErrEig(patchRotErrEig ~= 0))/pi*180) > 10) / sum(sum(patchRotErrEig ~= 0));

statistics = full([tau ; epsZinp ; epsZeig ; epsRinp ; epsReig ; outliersRinp ; outliersReig]);





