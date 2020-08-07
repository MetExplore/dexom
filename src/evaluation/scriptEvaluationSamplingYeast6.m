% Script for reproducing the evaluation using the DAG example
% https://www.biorxiv.org/content/10.1101/2020.07.17.208918v1

dexomInit;

% TESTED ONLY WITH IBM_CPLEX
% changeCobraSolver('ibm_cplex','all',0);

% CONFIGURE GLOBAL LP AND MILP OPTIONS
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit',30);
changeCobraSolverParams('MILP','logFile',0);
changeCobraSolverParams('LP','logFile',0);

load('gem_yeast6_06_reduced.mat');
maxSolutions = 1000;

% CONFIGURE GLOBAL ENUM AND METHOD OPTIONS
e.maxEnumTime = 30 * maxSolutions; % 30s max per solution on average
e.maxIterations = 1500;
e.maxSolutions = maxSolutions;
e.metricsUpdateFrequency = 0;
e.maxUniqueSolutions = maxSolutions;

nRxn = 60;

randData.RH = randsample(1:length(model.rxns), nRxn);
randData.RL = randsample(setdiff(1:length(model.rxns), randData.RH), nRxn);

m.distAnnealing = 0.995;
m.useRandomSeed = 1;
m.RHindex = randData.RH;
m.RLindex = randData.RL;

methods = {'dexom-default','dexom-rxnenum','dexom-icut','dexom-maxdist'}; 

results = {};
numRuns = 10;

basePath = fileparts(which('yeast_output_eval.txt'));
folderName = datestr(now,'dd-mm-yyyy');
outputPath = [basePath filesep folderName];

if ~exist(outputPath)
    mkdir(outputPath);
end


for j = 1:numRuns
    fprintf('run %d \n', j);
    for i = 1:length(methods)
        m.name = methods{i};
        o = setupMethodOptions(m.name, m, e);
        results{i, j} = sequentialNetworkEnumeration(model, o);
    end
    %save([outputPath filesep 'temporalResults'], 'results', '-v7.3');
end

outputResults = struct;
outputResults.results = results;
outputResults.methods = methods;
outputResults.model = model;
save([outputPath filesep 'results'], 'outputResults', '-v7.3');
