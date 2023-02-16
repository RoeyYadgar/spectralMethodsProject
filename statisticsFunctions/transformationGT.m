function [reflections,rotations] = transformationGT(patches,posGT)
%finds the ground truth tranformation of each patch
%inputs:
%   patches - N x 2 cell matrix of patches (first column) and their
%       corresponding star node (second column)
%   posGT - M x 2 ground truth position of all nodes
%outputs:
%   reflecitons - Nx1 vector of GT reflection of each patch
%   rotations - Nx1 vector of GT rotations of each patch

N = length(patches);
reflections = zeros(N,1);
rotations = zeros(N,1);
for i = 1:N
    patchGT = patches{i,1};
    ID = patchGT.Nodes.ID;
    patchGT.Nodes.Pos = posGT(ID,1) + 1i*posGT(ID,2);
    [ref,rot] = alignPatchesLSregis(patchGT,patches{i,1},ID);
    reflections(i) = ref;
    rotations(i) = conj(rot); 
end




end