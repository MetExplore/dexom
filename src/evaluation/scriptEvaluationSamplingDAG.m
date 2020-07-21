% Script for reproducing the evaluation using the DAG example
% https://www.biorxiv.org/content/10.1101/2020.07.17.208918v1

% Make sure that DEXOM is initialized
dexomInit;

% Set to 1 to use integer-cuts to perform full enumeration
fullEnumeration = 0;

basePath = fileparts(which('evaluation.md'));
basePath = [basePath filesep 'results' filesep 'dag'];
folderName = datestr(now,'dd-mm-yyyy');
outputPath = [basePath filesep folderName];

if fullEnumeration == 1
    outputPath = [basePath filesep folderName filesep 'full'];
end

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

logFile = [outputPath filesep 'output.log'];
diary(logFile);
diary ON;


% TESTED ONLY WITH IBM_CPLEX
% changeCobraSolver('ibm_cplex','all',0);

% CONFIGURE GLOBAL LP AND MILP OPTIONS
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit',30);
changeCobraSolverParams('MILP','logFile',0);
changeCobraSolverParams('LP','logFile',0);

methods = {'dexom-default','dexom-rxnenum','dexom-icut','dexom-maxdist'}; 


numLayers = 5;
numMetsPerLayer = 4;
model = dagNet(numLayers, numMetsPerLayer);
maxSolutions = 250;
numRuns = 30;

if fullEnumeration == 1
    methods = {'dexom-icut'};
    maxSolutions = numMetsPerLayer^numLayers;
    maxIterations = maxSolutions;
    maxUniqueSolutions = maxSolutions;
    numRuns = 1;
end

% CONFIGURE GLOBAL ENUM AND METHOD OPTIONS
e.maxEnumTime = 30 * maxSolutions; % 30s max per solution on average
e.maxIterations = maxSolutions * 3;
e.maxSolutions = maxSolutions;
e.metricsUpdateFrequency = 0;
e.maxUniqueSolutions = maxSolutions;

m.RHindex = [];
m.RLindex = 1:length(model.rxns);
m.useRandomSeed = 1;


results = {};
for j = 1:numRuns
    fprintf('run %d \n', j);
    for i = 1:length(methods)
        m.name = methods{i};
        o = setupMethodOptions(m.name, m, e);
        results{i, j} = sequentialNetworkEnumeration(model, o);
        [numSols, ~] = size(results{i,j}.solutions);
        maxSolutions = min(maxSolutions, numSols);
    end
end


outputResult = struct;
outputResult.results = results;
outputResult.methods = methods;
outputResult.model = model;

save([outputPath filesep 'results'], 'outputResult', '-v7.3');
exportResults(results, outputPath);

diary OFF;
