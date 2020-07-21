function classifications = evaluateModels(model, koFile, growthResults, aggregationMethod, koTol)
    koTable = readtable(koFile);
    % Iterate over conditions (i) and thresholds (j)
    for i = 1:size(growthResults, 1)
        for j = 1:size(growthResults, 2)
            % If there are multiple optimal networks, use the ensemble. By
            % default use the min function to reduce the networks (i.e., if
            % one network predicts a gene as essential, 
            growthRatios = growthResults{i, j};
            if size(growthRatios, 1) > 1
                switch aggregationMethod
                    case 'min'
                        growthRatios = min(growthRatios);
                    case 'max'
                        growthRatios = max(growthRatios);
                    case 'mean'
                        growthRatios = mean(growthRatios);
                    case 'first'
                        % Consider only the first sol (default solution)
                        growthRatios = growthRatios(1,:);
                    otherwise
                        error('Unrecognized aggregation method');
                end
            end
            classifications{i,j} = evaluateKO(model, koTable, growthRatios', koTol);
        end
    end
end