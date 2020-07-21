function [model, data] = mapYeastData(cobraModel, ...
    dataset, condition, quantiles, minGrowthFlux, constraintMode)
    
    model = configureYeastModel(cobraModel, constraintMode, minGrowthFlux);
    condIdx = find(strcmp(condition, dataset.conditions));
    transcriptomics = dataset.transcriptomics(:, condIdx);
    expression.gene = dataset.genes;
    expression.value = transcriptomics;
    
    % Use global thresholding to classify genes into highly expressed
    % (above upper threshold) or lowly expressed (below lower threshold)
    %lowerThreshold = quantile(dataset.transcriptomics(:), quantiles(1));
    %upperThreshold = quantile(dataset.transcriptomics(:), quantiles(2));
    lowerThreshold = quantile(transcriptomics, quantiles(1));
    upperThreshold = quantile(transcriptomics, quantiles(2));

    reactionExpression = mapExpressionToReactions(model, expression);
    reactionExpression(reactionExpression < 0) = nan; % unknown

    RH = reactionExpression >= upperThreshold;
    RL = reactionExpression < lowerThreshold;
    
    reactionLevels = reactionExpression;
    reactionLevels(RH) = 1;
    reactionLevels(RL) = -1;
    reactionLevels(~(RH|RL)) = 0;
    
    data.RHindex = find(RH);
    data.RLindex = find(RL);
    data.reactionExpression = reactionExpression;
    data.reactionLevels = reactionLevels;
    data.lowerThreshold = lowerThreshold;
    data.upperThreshold = upperThreshold;
end

