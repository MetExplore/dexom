%dexomInit;

% TESTED ONLY WITH IBM_CPLEX
changeCobraSolver('ibm_cplex','all',0);


%% LOG
path = fileparts(which('evaluation.md'));
folderName = datestr(now,'dd-mm-yyyy');
outputPath = [path filesep 'results' filesep 'cancer' filesep folderName];

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

globalOptions = struct;
globalOptions.solverName = 'ibm_cplex';

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
% flux for highly express reactions
globalOptions.method.epsilon = 1e-4;
globalOptions.method.runtime = globalOptions.solver.timeLimit;
globalOptions.method.printLevel = 0;
globalOptions.method.distAnnealing = 0.995;

% ENUM OPTIONS
globalOptions.enum.metricsUpdateFrequency = 60;
globalOptions.enum.maxIterations = 3300;
globalOptions.enum.maxEnumTime = 60 * 60 * 8;
% Optimal tolerance for accepting solutions during enumeration
%globalOptions.enum.optTol = globalOptions.solver.relMipGapTol * 1.5;
globalOptions.enum.optTol = 0.01;
globalOptions.enum.verbose = 1;
globalOptions.enum.maxTries = 20;
globalOptions.enum.maxUniqueSolutions = 3000;
globalOptions.enum.maxSolutions = 3300;

%globalOptions.enum.maxSolutions = 10;

% Configure the solver
changeCobraSolver(globalOptions.solverName,'all',0);
changeCobraSolverParams('LP','printLevel',0);
changeCobraSolverParams('LP','logFile',0);
changeCobraSolverParams('MILP','relMipGapTol',globalOptions.solver.relMipGapTol);
changeCobraSolverParams('MILP','printLevel',0);
changeCobraSolverParams('MILP','timeLimit',globalOptions.solver.timeLimit);
changeCobraSolverParams('MILP','logFile',0);

%% LOAD DATA

cell = 'A375';

load(strcat('gene_expr_u_', cell));
load(strcat('gene_id_u_', cell));
load(strcat('ID_FPKM_', cell));
load(strcat('model_u_', cell));


pseudocounts = log10(1 + gene_expr_u);
q1 = quantile(pseudocounts, 0.10);
q2 = quantile(pseudocounts, 0.90);

% Use min/max to parse data (OR GPR -> MAX, AND GPR -> MIN)
% Use 1 for highly expressed, -2 for lowly expressed. -1 is used
% for unknown.
expressionData.gene = gene_id_u;   
expressionData.value = pseudocounts;
    
fprintf('Mapping gene expression to reactions... ');
[mapping, ~, ~] = mapExpressionToReactions(model_u, expressionData);
fprintf('OK\n');

results.mapping = mapping;
o = setupMethodOptions('dexom-default/norand', globalOptions.method, globalOptions.enum);

o.method.RHindex = find(mapping > q2 );
o.method.RLindex = find(mapping >= 0 & mapping <= q1);
fprintf('RH %d, RL %d\n', length(o.method.RHindex), length(o.method.RLindex));

result = sequentialNetworkEnumeration(model_u, o);
        
acceptedFluxSolutions = result.contSolutions(result.accepted == 1, :);
active = abs(acceptedFluxSolutions) >= globalOptions.method.tol;
innactive = ~active;
acceptedFluxSolutions(active) = 1;
acceptedFluxSolutions(innactive) = 0;
uniqueDiscretized = unique(acceptedFluxSolutions, 'rows');
results.unique = uniqueDiscretized;
[enrichment,pathways] = pathwayEnrichmentDistribution(model_u, results);

diary OFF;