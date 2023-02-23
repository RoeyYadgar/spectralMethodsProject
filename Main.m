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
    title(['Reconstruction of the data points with the ASAP algorithm' newline 'Eta = ' num2str(100*noiseLevel(i)) ' %']);
    [aneScore,alignmentStats] = asapAlg(dataPoints,rho,noiseLevel(i));
    resultTable(i,'NoiseLevel') = {noiseLevel(i)};
    resultTable(i,'ANEScore') = {aneScore};
    resultTable(i,3:end) = alignmentStats';
end
