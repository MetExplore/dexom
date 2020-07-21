function options = setupMethodOptions(methodName, methodOptions, enumOptions)
    % Get default options and overwrite with user-specified options
    % enumOptions are optional
    if ~exist('enumOptions','var'), enumOptions = struct; end 
    
    % Generic, common options
    mo = updateAttributes(defaultMethodOptions(), methodOptions, 1);
    eo = updateAttributes(defaultEnumOptions(), enumOptions, 1);
    mo.name = methodName;
    
    % Split to get the initial part
    methodParts = strsplit(mo.name, '-');
    
    % Overwrite specific options depending on the method
    switch methodParts{1}
        case 'dexom'
            mo = dexomSetupOptions(mo);
            eo.enumGreedySkipReactions = 1;
            mo.useDefaultIterativeModelConstraintUpdate = 0;
            switch methodParts{2}
                case 'icut'
                    eo.enumStrategy = 'none';
                    mo.useSecondaryObjective = 0;
                    mo.excludeSolutionsByIntegerCuts = 1;
                case 'maxdist'
                    eo.enumStrategy = 'none';
                    mo.useSecondaryObjective = 1;
                    mo.secondaryObjectiveMode = 'maxdist';
                    mo.excludeSolutionsByIntegerCuts = 1;
                case 'mindist'
                    % Maximize overlap
                    eo.enumStrategy = 'none';
                    mo.secondaryObjectiveSense = -1;
                    mo.useSecondaryObjective = 1;
                    % Same as maxdist, but maximizing overlapping instead
                    % of minimizing it
                    mo.secondaryObjectiveMode = 'maxdist';
                    mo.excludeSolutionsByIntegerCuts = 1;
                case 'default'
                    % Start with default enum
                    eo.enumStrategy = 'random';
                    mo.useSecondaryObjective = 1;
                    mo.secondaryObjectiveMode = 'dexom';
                    mo.excludeSolutionsByIntegerCuts = 1;
                case 'default/norand'
                    eo.enumStrategy = 'default';
                    mo.useSecondaryObjective = 1;
                    mo.secondaryObjectiveMode = 'dexom';
                    mo.excludeSolutionsByIntegerCuts = 1;    
                case 'default/nocuts'
                    % Start with default enum
                    eo.enumStrategy = 'none';
                    mo.useSecondaryObjective = 1;
                    mo.secondaryObjectiveMode = 'dexom';
                    mo.excludeSolutionsByIntegerCuts = 0;
                case 'uniform'
                    eo.enumStrategy = 'none';
                    mo.useSecondaryObjective = 1;
                    mo.secondaryObjectiveMode = 'random';
                    mo.excludeSolutionsByIntegerCuts = 0;
                case 'rxnenum'
                    mo.useSecondaryObjective = 0;
                    mo.excludeSolutionsByIntegerCuts = 0; 
                    eo.enumGreedySkipReactions = 0;
                    eo.enumStrategy = 'random';
                    mo.useDefaultIterativeModelConstraintUpdate = 1;
                otherwise
                    error('Unknown DEXOM method');
            end
            
        case 'cobra'
            switch methodParts{2}
                case 'imat'
                    eo.enumGreedySkipReactions = 0;
                    eo.enumStrategy = 'default';
                    mo.useDefaultIterativeModelConstraintUpdate = 1;
                case 'imat/random'
                    eo.enumGreedySkipReactions = 0;
                    eo.enumStrategy = 'random';
                    mo.useDefaultIterativeModelConstraintUpdate = 1;
                case 'imat/greedyskip'
                    eo.enumGreedySkipReactions = 1;
                    eo.enumStrategy = 'default';
                    mo.useDefaultIterativeModelConstraintUpdate = 1;   
                otherwise
                    error('No available cobra method');
            end
        otherwise
            error('Unrecognized method');
    end
    
    options.enum = eo;
    options.method = mo;
end

