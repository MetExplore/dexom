function exportResults(results, outputPath)   
    [methods, repetitions] = size(results);
    
    maxSolutions = 0;
    for j = 1:repetitions
        for i = 1:methods
            [numSols, ~] = size(results{i,j}.solutions);
            maxSolutions = max(maxSolutions, numSols);
        end
    end
    
    meanDistances = zeros(methods * repetitions, maxSolutions);
    meanNN = zeros(methods * repetitions, maxSolutions);
    elapsedTime = zeros(methods * repetitions, maxSolutions);
    labels = [];
    solutionLabels = [];
    solutionLabels2 = [];
    allSolutions = [];
    allSolutionsFromFlux = [];
    idx = 1;
    for i = 1:methods
        methodAllSols = [];
        methodAllSolsFromFlux = [];
        for j = 1:repetitions
			meanDist = results{i, j}.metrics.distance.mean;
			meanDistNN = results{i, j}.metrics.distance.avgNearest;
			elapsed = results{i, j}.elapsedTime;
			vecMeanDist = ones(1, maxSolutions) * meanDist(end);
			vecMeanDistNN = ones(1, maxSolutions) * meanDistNN(end);
			vecTime = ones(1, maxSolutions) * elapsed(end);
			vecTime(1:numel(elapsed)) = elapsed;
			vecMeanDist(1:numel(meanDist)) = meanDist;
			vecMeanDistNN(1:numel(meanDistNN)) = meanDistNN;
            meanDistances(idx, :) = vecMeanDist;
            meanNN(idx, :) = vecMeanDistNN;
            elapsedTime(idx,:) = vecTime;
            labels = [labels string(results{i, j}.options.method.name)];
            methodSols = results{i, j}.intSolutions(results{i, j}.intAccepted == 1,:);
            methodSolsFromFlux = results{i, j}.solutions(results{i, j}.accepted == 1,:);
            methodAllSolsFromFlux = [methodAllSolsFromFlux; methodSolsFromFlux];
            methodAllSols = [methodAllSols; methodSols];
            idx = idx + 1;
        end
        un = abs(round(unique(methodAllSols, 'rows')));
        un2 = abs(round(unique(methodAllSolsFromFlux, 'rows')));
        solutionLabels = [solutionLabels repelem(string(results{i, j}.options.method.name), size(un, 1))];
        solutionLabels2 = [solutionLabels2 repelem(string(results{i, j}.options.method.name), size(un2, 1))];
        allSolutions = [allSolutions; un];
        allSolutionsFromFlux = [allSolutionsFromFlux; un2];
    end
    
    fid = fopen([outputPath filesep 'allSolutionLabels.csv'],'w');
    fprintf(fid,'%s\n',solutionLabels');
    fclose(fid);
    
    fid = fopen([outputPath filesep 'allSolutionLabelsFromFlux.csv'],'w');
    fprintf(fid,'%s\n',solutionLabels2');
    fclose(fid);
    
    fid = fopen([outputPath filesep 'dist_labels.csv'],'w');
    fprintf(fid,'%s\n',labels');
    fclose(fid);
    
    csvwrite([outputPath filesep 'allSolutionsIntVar.csv'], allSolutions);
    csvwrite([outputPath filesep 'allSolutionsFromFlux.csv'], allSolutionsFromFlux);
    csvwrite([outputPath filesep 'mean_nn_dist.csv'], meanNN);
    csvwrite([outputPath filesep 'mean_dist.csv'], meanDistances);

end