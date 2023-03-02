%% Globally rigid example
figure
nodes = zeros(5,2);
nodes(1,:) = [-1 -1]; nodes(2,:) = [-1,1] ; nodes(4,:) = [1,1]; nodes(5,:) = [1,-1];
c = zeros(5,5); c(:,3) = 1; c(3,:) = 1; c(3,3) = 0; c(1,2) = 1; c(2,1) = 1; c(4,5) = 1; c(5,4) = 1;
p = graph(c);
p2 = graph(ones(5,5)-eye(5));
Rot = [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
nodes2 = nodes; nodes2(4:5,:) = Rot*nodes(4:5,:);
subplot(2,2,1)
plotPatch(p,nodes);
axis('equal')
title('Non Globally Rigid Graph')
subplot(2,2,3)
plotPatch(p,nodes2);
axis('equal')
subplot(1,2,2)
plotPatch(p2,nodes2);
axis('equal')
title('Globally Rigid Graph')
%% Globally Rigid Components of Star Graph Example
nodes = [0 0; -1 0; -1 -1; 0 1; 1 1; 1 0; 1 -0.5; 0.5 -1];
v1 = [ones(1,7) [2 2 3 4 5 6 6 7]];
v2 = [2:8 [3 4 4 5 6 7 8 8]];
g = graph(v1,v2,1:length(v1),{'1','2','3','4','5','6','7','8'});
g.Nodes.Pos = nodes;
patches = starGraphRigidPatches(g);
figure 
subplot(2,1,1)
plotPatch(g)
title('Globally Rigid Components of Star Graph Example')
subplot(2,1,2)
hold on
for i = 1:length(patches)
    pos = patches{i}.Nodes.Pos;
    plotPatch(patches{i},pos + mean(pos))
end
%% Combinatorial Score Method Example
rng('default');
rng(2);
n = 3;
nodes1 = rand(n,2);
nodes2 = rand(n,2); nodes2(:,1) = -nodes2(:,1);
commonNodes = [0 0.2 ; 0 0.8];
nodes = [nodes1;commonNodes; nodes2];
distances = sqrt((nodes(:,1) - nodes(:,1)').^2 + (nodes(:,2) - nodes(:,2)').^2);
adj = (distances <= 0.5) - eye(length(distances));
g = graph(adj);
ind1 = 1:(n+2);
ind2 = (n+1):(2*n+2);
g1 =graph(adj(ind1,ind1));
g2 = graph(adj(ind2,ind2));
figure
subplot(2,1,1)
plotPatch(g,nodes)
hold on
plotPatch(g1,nodes(ind1,:))
plotPatch(g2,nodes(ind2,:))
legend('Cross Patch Edges','Edges in Patch 1','Edges in Patch 2')
title(['Combinatorial Score Method Example' newline 'Correct Reflection'])

subplot(2,1,2)
nodes(ind2,1) = -nodes(ind2,1);
plotPatch(g,nodes)
hold on
plotPatch(g1,nodes(ind1,:))
title('Wrong Reflection')
plotPatch(g2,nodes(ind2,:))