% Script for reproducing the evaluation using the DAG example
% https://www.biorxiv.org/content/10.1101/2020.07.17.208918v1

% Make sure that DEXOM is initialized
dexomInit;

changeCobraSolver('ibm_cplex','all',0);
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('MILP','relMipGapTol', 1e-6);
changeCobraSolverParams('MILP','intTol', 1e-8);
changeCobraSolverParams('MILP','feasTol', 1e-9);
changeCobraSolverParams('MILP','optTol', 1e-9);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit', 30);
changeCobraSolverParams('MILP','logFile',0);
changeCobraSolverParams('LP','logFile',0);

load('evaluation/data/yeast/gem_yeast6_06_reduced');

nRxn = 60;

enumOpts.maxEnumTime = 60 * 10;
enumOpts.optTol = 0.025;

randData.RH = randsample(1:length(model.rxns), nRxn);
randData.RL = randsample(setdiff(1:length(model.rxns), randData.RH), nRxn);

methodOptions.RHindex = randData.RH;
methodOptions.RLindex = randData.RL;

methodOptions.epsilon = 1e-4;
methodOptions.tol = 1e-8;
methodOptions.decreaseFitScoreAfterTries = 3;

methods = {'dexom-default';'dexom-rxnenum';'dexom-maxdist';'dexom-icut'};
       
results = {};
for i = 1:length(methods)
    methodOptions.name = methods{i};
    o = setupMethodOptions(methodOptions.name, methodOptions, enumOpts);
    results{i} = sequentialNetworkEnumeration(model, o);
end


figure;
title('MEAN');
hold on;
hs = [];
for i = 1:length(results)
    hs(i) = plot(results{i}.elapsedTime, results{i}.metrics.distance.mean);    
end
legend(hs, methods);
hold off;

figure;
title('AVG NEAREST');
hold on;
hs = [];
for i = 1:length(results)
    hs(i) = plot(results{i}.elapsedTime, results{i}.metrics.distance.avgNearest);
end
legend(hs, methods);
hold off;

% Unique solutions
uniqueSolutions = {};
uniqueSolFromFlux = {};
fullMatrix = [];
fullMatrixS = [];
labels = [];
cLabels = {};
for i = 1:length(results)
    uniqueSolutions{i} = results{i}.intSolutions(results{i}.intAccepted == 1, :);
    uniqueSolFromFlux{i} = results{i}.solutions(results{i}.accepted == 1, :);
    cLabels{i} = repmat(methods(i),size(uniqueSolutions{i},1),1);
    labels = [labels; cLabels{i}];
    fullMatrix = [fullMatrix; uniqueSolutions{i}];
    fullMatrixS = [fullMatrixS; uniqueSolFromFlux{i}];
    d = pdist(uniqueSolutions{i}, 'hamming');
    figure;
    title(methods(i));
    hold on;
    hist(d, 1000);
    hold off;
end


% PCA PLOT
% Compute also the medoid of all the techniques
distances = pdist(fullMatrix, 'hamming');
Z = squareform(distances);
totalDist = sum(Z,2);
meanDist = mean(Z,2);
[~,medoidID] = min(totalDist);
[~,mDistantID] = max(totalDist);
medoid = fullMatrix(medoidID,:);

%[loadings, scores] = pca(fullMatrix);
[loadings,scores,latent,tsquared,explained,mu] = pca(fullMatrix);

figure;
hold on;
gscatter(scores(medoidID,1), scores(medoidID, 2), labels(medoidID), 'k', '*', 20);
gscatter(scores(mDistantID,1), scores(mDistantID, 2), labels(mDistantID), 'k', 'x', 20);
gscatter(scores(:,1), scores(:,2), labels, [], 'o')
hold off;


