function grRatios = conditionKnockouts(model, results, createUnionNetwork, verbose)
    % Each row is a condition 
    for i = 1:size(results, 1)
        % Each column is a threshold
        for j = 1:size(results, 2)
            r = results{i,j};
            fprintf('Running condition %s with thresholds [%f %f]...\n', r.options.condition, r.options.thresholds(1), r.options.thresholds(2));
            networks = r.uniqueOptimalNetworks;
            m = model;
            m.lb = r.model_lb;
            m.ub = r.model_ub;
            m.lb(find(strcmp(m.rxnNames, 'growth'))) = 0;
            grRatios{i,j} = multiNetKnockouts(m, networks, createUnionNetwork, verbose);
        end
    end
end