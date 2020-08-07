function results = runExperiments(model, dataset, conditions, methods, thresholds, ...
    mediums, biomasses, globalOptions, storeFolder)
    global CBT_LP_SOLVER;
    
    storeOutput = 0;
    if exist('storeFolder','var')
        if exist(storeFolder, 'dir')
            storeOutput = 1;
        else
            fprintf('Output folder does not exist, ignoring folder');
        end
    else
        warning('Output is not going to be saved!');
    end
    % Generate the experimental matrix design
    % Condition, Medium, Biomass, Threshold, Method
    p = cartesianProd(1:size(conditions,1), 1:size(mediums,1), ...
        1:size(biomasses,1), 1:size(thresholds,1), 1:size(methods,1));
    
    tableData = {};
    
    % Filename to store each result
    basename = datestr(now,'dd-mm-yyyy-HH-MM');
    basename = strcat('yeastResult_', basename);

    for i = 1:size(p,1)
        % For each setting, reconstruct networks and evaluate the
        % performance on gene essentiality prediction
        
        expSettings = p(i,:);
        info.condition = conditions{expSettings(1)};
        info.medium = mediums{expSettings(2)};
        info.minBioFlux = biomasses(expSettings(3));
        info.thresholds = thresholds(expSettings(4),:);
        info.method = methods{expSettings(5)};
        printInfo(info);
        
        switch info.medium
            case 'full'
                mediumNumber = 0;
            case 'minimal'
                mediumNumber = 1;
            case 'anaerobic'
                mediumNumber = 2;
            case 'default'
                mediumNumber = -1;
        end
        
        % Prepare the data for this setting
        [m, data] = mapYeastData(model, dataset, info.condition, ...
            info.thresholds, info.minBioFlux, mediumNumber);
        
        % Configure method (wrap this to do extra checks and skip configs)
        o = setupMethodOptions(info.method, globalOptions.method, globalOptions.enum);
        o.method.RHindex = data.RHindex;
        o.method.RLindex = data.RLindex;
        result = sequentialNetworkEnumeration(m, o);
        
        % TODO: Save algorithm information
        if isfield(result, 'skipped') && result.skipped == 1
            fprintf('Skipped condition\n');
            continue;
        end
        
        % Reset growth bounds before KO testing
        biomassId = find(strcmp(m.rxnNames, 'growth'));
        m.lb(biomassId) = 0;
        
        % Recover only the accepted, unique solutions (above optimal tol)
        solutions = unique(full(result.solutions(result.accepted == 1,:)), 'rows');
        numSolutions = size(solutions,1);
        countRxn = sum(solutions, 2);
        
        % Evaluate KOs in each network, add row per solution, if more than
        % 1 sol, show also aggregations
        % TODO: Move evaluation of the networks outside the method
        fprintf('\t- Number of solutions: %d\n', numSolutions);
        fprintf('\t- Number of reactions: mean %f, sd %f\n', mean(countRxn), std(countRxn));

        if ~isfield(dataset, 'koDataset') || isempty(dataset.koDataset)
            fprintf('No dataset for KO evaluation provided, skipping gene deletion\n');
            continue;
        end
        
        % SIMULATE KOs
        result.evalNetsKO = cell(1,numSolutions);
        result.growthRatios = zeros(numSolutions, length(m.genes));
        bestTPR = 0;
        worstFPR = 0;
        
        koD = dataset.koDataset;
        koTol = globalOptions.koTol;
        cgRatios = cell(1, numSolutions);
        cScores = cell(1, numSolutions);
        
        environment = getEnvironment();
        parfor j = 1:numSolutions
            restoreEnvironment(environment);
            changeCobraSolver(CBT_LP_SOLVER, 'LP', 0, -1); 
            s = solutions(j,:);
            fprintf('Simulating WT vs KO growth on network %04d (%d rxn) ... ', j, sum(s));
            grRatios = simulateKO(m, s); 
            score = evaluateKO(m, koD, grRatios, koTol);
            fprintf('Done (%.4f TPR, %.4f FPR, %.4f ACC, %.4f MCC)\n', ...
                score.TPR, score.FPR, score.accuracy, score.mcc);
            cgRatios{j} = grRatios;
            cScores{j} = score;
        end
        
        for j = 1:numSolutions
            tableData = addResult(tableData, info, 'single', cScores{j});
            result.evalNetsKO{j} = cScores{j};
            result.growthRatios(j,:) = cgRatios{j};
        end
        
        %  EVALUATRE KOs
        result.evalUnionNetKO = [];
        result.evalOrKO = [];
        result.evalAndKO = [];
        
        if numSolutions > 1
            % Union Network
            unionNet = max(solutions,[],1);
            grRatiosUnion = simulateKO(model, unionNet);
            result.evalUnionNetKO = evaluateKO(model, dataset.koDataset, grRatiosUnion, globalOptions.koTol);
            tableData = addResult(tableData, info, 'union', result.evalUnionNetKO, 1);
            % Or Ensemble
            minGrowthPrediction = min(result.growthRatios,[],1);
            result.evalOrKO = evaluateKO(model, dataset.koDataset, minGrowthPrediction', globalOptions.koTol);
            tableData = addResult(tableData, info, 'or', result.evalOrKO, 1);
            % And Ensemble
            maxGrowthPrediction = max(result.growthRatios,[],1);
            result.evalAndKO = evaluateKO(model, dataset.koDataset, maxGrowthPrediction', globalOptions.koTol);
            tableData = addResult(tableData, info, 'and', result.evalAndKO, 1);
        end
        
        store.result = result;
        store.info = info;
        filename = strcat(basename, sprintf('_method-%s-th-%.2f-%.2f-ex-%d.mat', info.method, info.thresholds(1), info.thresholds(2), i));
        if storeOutput
            save([storeFolder filesep filename], 'store');
        end
        
        if storeOutput
            save([storeFolder filesep strcat('table-',basename,'.mat')],'tableData','-v7.3');
        end
    end
    
    if isfield(dataset, 'koDataset') && ~isempty(dataset.koDataset)
        tableData = cell2table(tableData,'VariableNames', ...
            {'Condition', 'Medium', 'MinBiomassFlux', 'Thresholds', ...
            'Method', 'Mode', 'MCC', 'F1', 'Accuracy', ...
            'TPR', 'FPR', 'TP', 'FP', 'TN', 'FN'});

        results.evaluation = tableData;
        if storeOutput
            save([storeFolder filesep strcat('table-',basename,'.mat')],'tableData','-v7.3');
        end
    end

end

function printInfo(info)
    fprintf('> Condition = %s, Method = %s, Threshold = [%.2f, %.2f], Medium = %s, Min. Biomass flux = %.4f\n', ...
           info.condition, info.method, info.thresholds, info.medium, info.minBioFlux);
end

function results = addResult(results, info, mode, score, showInfo)
    %printInfo(info);
    %fprintf('\t > Evaluation result (mode %s): MCC = %f, F1 = %f, Accuracy = %f, TPR = %f, FPR = %f\n', ...
    %    mode, score.mcc, score.f1, score.accuracy, score.sensitivity, 1 - score.specificity);
    if nargin < 5
        showInfo = 0;
    end
    
    if showInfo == 1
        fprintf('Evaluation of the %s-Network: %.4f TPR, %.4f FPR, %.4f ACC, %.4f MCC\n', ...
            upper(mode), score.TPR, score.FPR, score.accuracy, score.mcc);
    end
    
    results(end+1,:) = {info.condition, info.medium, info.minBioFlux, ...
        info.thresholds, info.method, mode, score.mcc, score.f1, ...
        score.accuracy, score.sensitivity, 1 - score.specificity, ...
        score.N_TP, score.N_FP, score.N_TN, score.N_FN};
end