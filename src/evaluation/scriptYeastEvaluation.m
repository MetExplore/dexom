%% LOG
path = fileparts(which('evaluation.md'));
folderName = datestr(now,'dd-mm-yyyy');
outputPath = [path filesep 'results' filesep 'yeast' filesep 'essential_genes' filesep folderName];

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

basename = datestr(now,'dd-mm-yyyy-HH-MM');
logFile = [outputPath filesep strcat(basename,'.log')];
resultFile = [outputPath filesep strcat(basename,'.mat')];
imageFile = [outputPath filesep strcat(basename,'.png')];

diary(logFile);
diary ON;

%% Parametes and configuration to run the reconstruction
studyDataFile = [path filesep 'data' filesep 'yeast' filesep 'rintalaYeast2009'];
modelFile = [path filesep 'data' filesep 'yeast' filesep 'gem_yeast6_06_reduced'];
%modelFile = [path filesep 'data' filesep 'yeast' filesep 'gem_yeast8_35_reduced'];
essentialDatasetFile = [path filesep 'data' filesep 'yeast' filesep 'dataset_yeast6_v1.0.7.csv'];

globalOptions = struct;
globalOptions.solverName = 'ibm_cplex';
% Simulated biomass flux ratio KO/WT < koTol is considered essential
globalOptions.koTol = 0.01;

% SPECIFIC SOLVER OPTIONS
globalOptions.solver.relMipGapTol = 1e-4;
%globalOptions.solver.relMipGapTol = 0.01;
% 5 min time limit per problem for the solver (if supported)
globalOptions.solver.timeLimit = 300; 

globalOptions.solver.optTol = 1e-6;

% METHOD OPTIONS
% Tolerance on flux values to include the reaction in the model. Once a
% solution is found, some reactions will carry a very small amount of flux
% below epsilon (so they are not considered highly expressed) but still
% needed to have a consistent steady state solution.
globalOptions.method.tol = 1e-6;
% flux for highly express reactions (iMAT, EXOM)
globalOptions.method.epsilon = 1e-4;
globalOptions.method.runtime = globalOptions.solver.timeLimit;
globalOptions.method.printLevel = 0;
globalOptions.method.distAnnealing = 0.995;

% ENUM OPTIONS
globalOptions.enum.metricsUpdateFrequency = 1e-8;
globalOptions.enum.maxIterations = 2500;
globalOptions.enum.maxEnumTime = 60 * 60 * 8;
% Optimal tolerance for accepting solutions during enumeration
%globalOptions.enum.optTol = globalOptions.solver.relMipGapTol * 1.5;
globalOptions.enum.optTol = 0.01;
globalOptions.enum.verbose = 1;
globalOptions.enum.maxTries = 100;
globalOptions.enum.maxUniqueSolutions = 2000;
globalOptions.enum.maxSolutions = 2500;

% Configure the solver
changeCobraSolver(globalOptions.solverName,'all',0);
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('LP','logFile',0);
changeCobraSolverParams('MILP','relMipGapTol',globalOptions.solver.relMipGapTol);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit',globalOptions.solver.timeLimit);
changeCobraSolverParams('MILP','logFile',0);


%% DATA LOAD

load(studyDataFile);
load(modelFile);

dataset.koDataset = readtable(essentialDatasetFile);

conditions = {'O2_20.9'};
          

thresholds = [0.10 0.90; 0.10 0.85; 0.10 0.80; 0.10 0.75;
              0.15 0.90; 0.15 0.85; 0.15 0.80; 0.15 0.75;
              0.20 0.90; 0.20 0.85; 0.20 0.80; 0.20 0.75;
              0.25 0.90; 0.25 0.85; 0.25 0.80; 0.25 0.75]; 
          
methods = {'dexom-default';'dexom-rxnenum';'dexom-maxdist';'dexom-icut'};
               
mediums = {'full'};
biomasses = 1e-4;

% TODO: Split method into generation of metabolic networks and evaluation
% with the gene essentiality dataset
results = runExperiments(model, dataset, conditions, ...
                         methods, thresholds, mediums, ...
                         biomasses, globalOptions, outputPath);
   
diary OFF;

%% PLOTS

hold on;
plots = [];
for i = 1:length(methods)
    r = results.evaluation(strcmp(results.evaluation.Method, methods{i}),:);
    plots(i) = scatter(r.FPR, r.TPR);
end
legend(plots, methods);
hold off;

saveas(gcf, imageFile)

%% STORE RESULTS
% save(resultFile, 'results', '-v7.3');
