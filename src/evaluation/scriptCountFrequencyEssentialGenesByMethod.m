% Path to the yeast6-evaluation-da0995-default.zip file
% Can be downloaded at https://doi.org/10.5281/zenodo.4279170
path = 'C:/work/dexom-results/yeast-essential-genes/yeast6-evaluation-da0995-default.zip';

t_diversity_enum = countFrequentEssentialGenesFromExportedFiles(path, '*dexom-default*');
t_reaction_enum = countFrequentEssentialGenesFromExportedFiles(path, '*dexom-rxnenum*');
t_maxdist_enum = countFrequentEssentialGenesFromExportedFiles(path, '*dexom-maxdist*');
t_icut_enum = countFrequentEssentialGenesFromExportedFiles(path, '*dexom-icut*');

% Load information about essentiality of each gene. Make sure that the
% dataset_yeast6_v.1.0.7.csv is in the path (it should be after calling
% dexomInit)
dataset = readtable('dataset_yeast6_v1.0.7.csv');

% Match names
idx = zeros(numel(t_diversity_enum.Var1),1);
for i = 1:numel(t_diversity_enum.Var1)
    idx(i) = find(strcmp(t_diversity_enum.Var1(i), dataset.GENE_ID));
end

essential = dataset(idx,:).LABEL_AEROBIC == 1;

tableFrequencies = table(t_diversity_enum.Var1, essential, t_diversity_enum.Var2, ...
    t_reaction_enum.Var2, t_maxdist_enum.Var2, t_icut_enum.Var2, ...
    'VariableNames',{'Gene','Essential','DiversityEnum','ReactionEnum','MaxdistEnum', 'IcutEnum'});