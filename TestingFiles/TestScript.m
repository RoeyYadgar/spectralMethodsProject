%Test script

v1 = [ones(1,7) [2 2 3 4 5 6 6 7 9]];
v2 = [2:8 [3 4 4 5 6 7 8 8 2]];
g = graph(v1,v2,1:length(v1),{'1','2','3','4','5','6','7','8','9'});
p = plot(g);
p.EdgeCData = biconncomp(g);
patches = starGraphRigidPatches(g);

%%
% load('USdistance.mat',['W']);
% rho = 0.32;
% W(W > rho) = 0;
% G = graph(W);

load('uscities.mat');
dataPoints = uscities(1:3,:)';
rho = 0.032;
G = generateGraphDiscModel(dataPoints,rho);
W = adjacency(G,'weighted');
patches = splitGraphToGloballiyRigidsComps(G);
%%
i = 2;
edges = table2array(G.Edges);
nodes1 = edges(find(edges(:,1) == i),2);
nodes2 = edges(find(edges(:,2) == i),1);
neighboringNodes = union(nodes1,nodes2);
nodeNames = cell(length(neighboringNodes)+1,1);
nodeNames{1} = num2str(i);
for j = 1:length(neighboringNodes)
    nodeNames{j+1} = num2str(neighboringNodes(j));
end
subStarGraph = subgraph(G,[i ;neighboringNodes]);
subStarGraph.Nodes = nodeNames;
patches = starGraphRigidPatches(subStarGraph);
p = patches{1};
x = localizePatch(p);
plot(p,'XData',x(:,1),'YData',x(:,2),'ZData',zeros(length(x),1));
figure
id = [i ; neighboringNodes];
plot(p,'XData',dataPoints(id,1),'YData',dataPoints(id,2),'ZData',dataPoints(id,2));
%%
N = 1907;
t = 2*(randi(2,N,1)-1.5);
gamma = sparse(diag(t)*B*diag(t.^-1));
T = sign(rand(N,N)-0.2);
z = gamma.*T;
Z = sparse(diag(sum(abs(z),2).^(-1)))*z;
[V,D] = eigs(Z,5)
plot(V(:,1))
%%
Y = (dataPoints(:,1)-dataPoints(:,1)').^2 + (dataPoints(:,2)-dataPoints(:,2)').^2 +(dataPoints(:,3)-dataPoints(:,3)').^2;
Y = sqrt(Y);
x = cmdscale(Y,2);
%%
N = length(patches);
reflections = zeros(N,1);
rotations = zeros(N,1);
for i = 1:N
    realPosPatch = patches{i,1};
%     realPosPatch.Nodes.Pos = 
    realPos = x(cellfun(@str2num,patches{i,1}.Nodes.Name),:);
    %[V,~] = eigs(realPos'*realPos,2);
    %realPosProj = realPos*V;
    realPosPatch.Nodes.Pos = realPos;
    [reflections(i),rotations(i)] = alignPatchesLSregis(realPosPatch,patches{i,1},realPosPatch.Nodes.Name);
end
%%
B = sparse(A >= 2);
pRef = diag(reflections)*B*diag(reflections);
pRot = diag(rotations)*B*diag(rotations.^-1);
refErr = pRef - patchReflection;
rotErr = pRot - patchRotation;
plot(sum(abs(refErr/2))./sum(abs(B)))
%%
subplot(2,1,1)
plotPatch(patches{i,1})
axis('equal')
subplot(2,1,2)
plotPatch(patches{i,1},x(cellfun(@str2num,patches{i,1}.Nodes.Name),:));
axis('equal')
%%
%Q1: what is the prescribed size in page 25?
%Q1.5: how to split a patch that is too big?
%Q2: unit measurements of rho and database
%Q3: use of matlab functions (mds, biconncomp)
%Q4: why cant you use the Combinatorial Score method to compute tranlation
%between patches?
%Q5: why cant you add every pair of nodes in each patch to the LS