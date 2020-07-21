% Load all the optimal networks enumerated by all methods
% for the threshold 0.25 0.90
load('gem_yeast6_06_reduced.mat');
[s,m] = combineMetabolicNetworksFromFile('yeast6-evaluation-nogit.zip', ...
    '07-06-2020/*-th-0.25-0.90*.mat');
% Get only the variable reactions (discard reactions that are always
% present or absent)
counts = sum(s, 1);
variableRxn = counts > 0 & counts < max(counts);
sv = s(:, variableRxn);
rxnNames = model.rxnNames(variableRxn);

% Export names of the variable reactions
path = [fileparts(which('dexom-evaluation.csv')) filesep 'networks'];
fpath = [path filesep 'reactionNames.txt'];
filePh = fopen('reactionNames.txt','w');
fprintf(filePh,'%s\n',rxnNames{:});
fclose(filePh);

% Export matrix
csvwrite([path filesep 'networksVarRxn.csv'], sv);

% Export method labels
filePh = fopen([path filesep 'methods.txt'],'w');
fprintf(filePh,'%s\n',m{:});
fclose(filePh);

