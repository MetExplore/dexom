% Load the Yeast 6 model used for the evaluation
load('gem_yeast6_06_reduced.mat');

% Load the dataset with the essential genes
dataset = 'dataset_yeast6_v1.0.7.csv';

% For each method, load the result files and aggregate the essential genes
fprintf('Evaluating DEXOM...\n');
dexomME = evaluateMetaEnsemble(model, 'yeast6-evaluation-nogit.zip', ...
    '07-06-2020/*-dexom-default*.mat', dataset);

fprintf('Evaluating Rxn-enum...\n');
rxnME = evaluateMetaEnsemble(model, 'yeast6-evaluation-nogit.zip', ...
    '07-06-2020/*-dexom-rxnenum*.mat', dataset);

fprintf('Evaluating Maxdist...\n');
maxdistME = evaluateMetaEnsemble(model, 'yeast6-evaluation-nogit.zip', ...
    '07-06-2020/*-dexom-maxdist*.mat', dataset);

fprintf('Evaluating integer-cut...\n');
icutME = evaluateMetaEnsemble(model, 'yeast6-evaluation-nogit.zip', ...
    '07-06-2020/*-dexom-icut*.mat', dataset);