function [reflection,rotation,unreflectedViolations,reflectedViolations] = alignPatchesCombScore(patch1,patch2,commonNodes,graphDistance,rho)
%computes rigid tranformation to align 2 patches that have 2 common
%nodes using Combinatorical Score Methods
%inputs:
%   patch1 - first patch
%   patch2 - second patch
%   commonNodes - array of common nodes between the two patches
%   graphDistance - distance matrix of the original graph
%   rho - distance threshold of the Disc Graph Model
%outputs:
%   reflection - reflection of the rigid tranformation (1 or -1)
%   rotation - rotation of the rigid transformation (e^(i*theta))

patch1Pos = patch1.Nodes.Pos(:,1) + 1i *  patch1.Nodes.Pos(:,2);
patch2Pos = patch2.Nodes.Pos(:,1) + 1i *  patch2.Nodes.Pos(:,2);
patch2ReflectedPos = conj(patch2Pos);

node1IndPatch1 = find(strcmp(commonNodes{1},patch1.Nodes.Name));
node1IndPatch2 = find(strcmp(commonNodes{1},patch2.Nodes.Name));
node2IndPatch1 = find(strcmp(commonNodes{2},patch1.Nodes.Name));
node2IndPatch2 = find(strcmp(commonNodes{2},patch2.Nodes.Name));



%Shift the patches such that node1 is at the origin in both patches
patch1Pos = patch1Pos - patch1Pos(node1IndPatch1);
patch2Pos = patch2Pos - patch2Pos(node1IndPatch2);
patch2ReflectedPos = patch2ReflectedPos - patch2ReflectedPos(node1IndPatch2);

%rotation to align the second node on the second patch to the first patch 
rotation = (patch1Pos(node2IndPatch1)/abs(patch1Pos(node2IndPatch1)))...
    *(patch2Pos(node2IndPatch2)/abs(patch2Pos(node2IndPatch2)))^-1;
rotationReflected = (patch1Pos(node2IndPatch1)/abs(patch1Pos(node2IndPatch1)))...
    *(patch2ReflectedPos(node2IndPatch2)/abs(patch2ReflectedPos(node2IndPatch2)))^-1;

patch2Pos = patch2Pos * rotation;
patch2ReflectedPos = patch2ReflectedPos * rotationReflected;

if(isnan(rotation) || isnan(rotationReflected))
    rotation = 0;
    reflection = 0;
    unreflectedViolations = 0;
    reflectedViolations = 0;
    return
end


unreflectedViolations = 0;
reflectedViolations = 0;
for i = 1:size(patch1Pos)
    for j = 1:size(patch2Pos)
        if(i ~= node1IndPatch1 && i ~= node2IndPatch1 && j~= node1IndPatch2 && j ~= node2IndPatch2)%only check for viloations when the nodes are neither of the common nodes
            edgeExists = graphDistance(str2double(patch1.Nodes.Name{i}),str2double(patch2.Nodes.Name{j})) > 0;
            nodesDistance = abs(patch1Pos(i) - patch2Pos(j));
            nodesDistanceReflected = abs(patch1Pos(i) - patch2ReflectedPos(j));
            
            unreflectedViolations = unreflectedViolations +     (edgeExists ~= (nodesDistance < rho));
            reflectedViolations = reflectedViolations + (edgeExists ~= (nodesDistanceReflected < rho));
            
            
        end
    end
end

if(unreflectedViolations == reflectedViolations)
    reflection = 0;
    rotation = 0;
elseif(unreflectedViolations < reflectedViolations)
    reflection = 1;
else
    reflection = -1;
    rotation = rotationReflected;
end

% plotPatch(patch1,patch1Pos);
% hold on
% plotPatch(patch2,patch2Pos);
% plotPatch(patch2,patch2ReflectedPos);

end
