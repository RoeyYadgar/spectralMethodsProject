load('uscities.mat');
dataPoints = uscities(1:3,:)';

%Data is on a sphere find its center and radius by performing iterations of linear approximation of the equation and solve with leastsquare
center = [0 0 0];
for i = 1:50
    dist = vecnorm(center-dataPoints,2,2);
    A = [(center-dataPoints)./dist -ones(length(dataPoints),1)];
    x = A\(-dist);
    center = center + x(1:3)';
end
R = x(4);

%Calculte arclength between points on the sphere
dataPoints = dataPoints - center;
dotProd = dataPoints*dataPoints';
dotProdNormalized = dotProd/(R^2);
dotProdNormalized(dotProdNormalized > 1) = 1;
theta = real(acos(dotProdNormalized));
W = R*theta;
W = W-diag(diag(W)); %set the diagonal to 0 (it's not exactly 0 due to numerial errors)

%Embed the data into a 2 dimensional space (so it can be plotted and compared to the estimated position)
dataPoints = cmdscale(W,2);

%%
rho = 0.032;
noiseLevel = 0:0.1:0.5;
resultTable = table('Size',[length(noiseLevel) 9],'VariableTypes',{'double','double','double','double','double','double','double','double','double'}...
    ,'VariableNames',{'NoiseLevel','ANEScore','Tau','EpsZInp','EpsZEig','EpsRInp','EpsREig','OutliersRinp','OutliersReig'});

figure
for i = 1:length(noiseLevel)
    subplot(2,3,i)
    [aneScore,alignmentStats] = asapAlg(dataPoints,rho,noiseLevel(i));
    resultTable(i,'NoiseLevel') = {noiseLevel(i)};
    resultTable(i,'ANEScore') = {aneScore};
    resultTable(i,3:end) = alignmentStats';
    if(i == 2)
        titletext = ['Reconstruction of the data points with the ASAP algorithm' newline 'Eta = ' num2str(100*noiseLevel(i)) ' %'];
    else
        titletext = ['Eta = ' num2str(100*noiseLevel(i)) ' %'];
    end
    title(titletext);
    drawnow
end
figure
uitable('Data',resultTable{:,:},'ColumnName',resultTable.Properties.VariableNames,...
    'RowName',resultTable.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 0.93],'fontSize',10);
title('Results of the ASAP algorithm')

%% Generating Example Figures used in the report
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
%%
nodes = [0 0; -1 0; -1 -1; 0 1; 1 1; 1 0; 1 -0.5; 0.5 -1];
v1 = [ones(1,7) [2 2 3 4 5 6 6 7]];
v2 = [2:8 [3 4 4 5 6 7 8 8]];
g = graph(v1,v2,1:length(v1),{'1','2','3','4','5','6','7','8'});
g.Nodes.Pos = nodes;
patches = starGraphRigidPatches(g);
figure 
subplot(2,1,1)
plotPatch(g)
title('Globally Rigid Components of Star Graph')
subplot(2,1,2)
hold on
for i = 1:length(patches)
    pos = patches{i}.Nodes.Pos;
    plotPatch(patches{i},pos + mean(pos))
end