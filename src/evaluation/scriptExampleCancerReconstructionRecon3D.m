%dexomInit;

%% Parametes and configuration to run the reconstruction
globalOptions = struct;

% SPECIFIC SOLVER OPTIONS
globalOptions.solver.relMipGapTol = 1e-4;
% 5 min time limit per problem for the solver (if supported)
globalOptions.solver.timeLimit = 300; 
globalOptions.solver.optTol = 1e-6;

% METHOD OPTIONS
% Tolerance on flux values to include the reaction in the model. Once a
% solution is found, some reactions will carry a very small amount of flux
% below epsilon (so they are not considered highly expressed) but still
% needed to have a consistent steady state solution.
globalOptions.method.tol = 1e-6;
globalOptions.method.runtime = globalOptions.solver.timeLimit;
globalOptions.method.printLevel = 0;
globalOptions.method.distAnnealing = 0.995;

% ENUM OPTIONS
globalOptions.enum.metricsUpdateFrequency = 1e-6;
globalOptions.enum.maxIterations = 20;
globalOptions.enum.maxSolutions = 20;
globalOptions.enum.maxEnumTime = 600;
globalOptions.enum.optTol = 0.01;
globalOptions.enum.verbose = 1;
globalOptions.enum.maxTries = 3;
globalOptions.enum.maxUniqueSolutions = 10;

% Configure the solver
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('LP','logFile',0);
changeCobraSolverParams('MILP','relMipGapTol',globalOptions.solver.relMipGapTol);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit',globalOptions.solver.timeLimit);
changeCobraSolverParams('MILP','logFile',0);

%% LOAD DATA

% This part requires the module 'modules/evaluation' included in the
% repository. This module is automatically added to the path if present
% when running the initialization script dexomInit.m.
cell = 'A375';
% Recon3D model is included in the evaluation repository, downloaded from:
% https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D_301.zip
load('Recon3DModel_301.mat');
model = Recon3DModel;

% Data from "A systematic evaluation of methods for tailoring genome-scale
% metabolic models" https://doi.org/10.1016/j.cels.2017.01.010
load(strcat('gene_expr_u_', cell));
load(strcat('gene_id_u_', cell));
load(strcat('ID_FPKM_', cell));

pseudocounts = log10(1 + gene_expr_u);
q1 = quantile(pseudocounts, 0.20);
q2 = quantile(pseudocounts, 0.80);

expressionData.gene = strcat(gene_id_u, '.1');   
expressionData.value = pseudocounts;
fprintf('Mapping gene expression to reactions... ');
[mapping, ~, ~] = mapExpressionToReactions(model, expressionData);
results.mapping = mapping;
fprintf('OK\n');

% Use any of the enumeration methods to get different reconstructions
% ('dexom-diversity', 'dexom-icut', 'dexom-reactionenum', 'dexom-maxdist')
o = setupMethodOptions('dexom-icut', globalOptions.method, globalOptions.enum);
o.method.RHindex = find(mapping > q2 );
o.method.RLindex = find(mapping >= 0 & mapping <= q1);
fprintf('RH %d, RL %d\n', length(o.method.RHindex), length(o.method.RLindex));

result = sequentialNetworkEnumeration(model, o);
solutions = getUniqueAcceptedSolutions(result);
