function tableResults = countFrequentEssentialGenesFromExportedFiles(zipFile, pattern)
    % Extract how frequent each gene was predicted as essential by each
    % method. It loads the raw results provided in a zip file
    % (yeast6-evaluation-da0995-default.zip) available at 
    % https://doi.org/10.5281/zenodo.4279170
    if ~exist('pattern', 'var') || isempty(pattern)
       pattern = '*dexom-default*';
   end 
    tmpFolder = tempname();
    mkdir(tmpFolder);
    unzip(zipFile, tmpFolder);
    % Load the Yeast 6.06 model
    load([tmpFolder filesep 'gem_yeast6_06_reduced.mat']);
    files = dir(fullfile(tmpFolder, strcat('**\', pattern)));
    z = zeros(1,900);
    t = 0;
    for i=1:numel(files)
        % Load stored data
        load([files(i).folder filesep files(i).name]);
        for j=1:numel(store.result.evalNetsKO)
            r = store.result.evalNetsKO{j};
            z = z + full(r.essentialGeneFromModel);
            t = t+1;
        end
    end
    freqPredictions = (z/t)';
    tableResults = table(model.genes, freqPredictions);
end