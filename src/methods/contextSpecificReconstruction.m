function solution = contextSpecificReconstruction(model, methodOptions, ctx)    
    % methodOptions.name should start by 'cobra-' (e.g. cobra-imat)
    % Common options:
    % - tol: lowest non-zero flux value
    % - epsilon: min flux for highly expressed reactions (iMAT)
    % - RHindex: indexes of the highly expressed (core) reactions (fastcore, iMAT, ...)
    % - RLindex: indexes of the lowly expressed reactions (ignored by other
    % methods)
    % - methodPrintLevel: value passed as printLevel to the method
    % - storeProblemDefinition: add the MILP or LP matrices to the result
    % output structure (.internal.problemDef, .internal.problemType)
    % - storeSolverOutput: add the output of solving the LP or MILP problem
    % (.internal.solverOutput)
    
    % All methods should return at least a boolRxnSolution vector with 
    % the reactions that carry non zero flux and the status of the solution
    % (valid, invalid)
    o = methodOptions;
    methodParts = strsplit(methodOptions.name, '-');
    
    switch methodParts{1}
        case 'cobra'
            switch methodParts{2}
                case {'imat','imat/greedyskip','imat/random'}
                    solution = cobraImat(model, o, ctx);
                    solution.activity = solution.boolRxnSolution;
                    solution.activity(solution.fluxSolution < 0 & solution.boolRxnSolution == 1) = -1;
                % TODO: add wrappers to other methods
                otherwise
                    error('Unsupported cobra method');
            end
        case 'dexom'
            solution = statefulDexom(model, o, ctx);
            solution.activity = solution.boolRxnSolution;
            solution.activity(solution.fluxSolution < 0 & solution.boolRxnSolution == 1) = -1;
        otherwise
            error('Unrecognized method type');                
    end
    
    % Depending on the objective function of the algorithm, the normalized
    % score would be different. Methods that only maximize and minimize RH
    % and RL sets should be evaluated using the coverage of both sets
    sizeRH = length(o.RHindex);
    sizeRL = length(o.RLindex);
    nRH = sum(solution.boolRxnSolution(o.RHindex) == 1);
    nRL = sum(solution.boolRxnSolution(o.RLindex) == 0);
    
    solution.normalizedScore = (nRH + nRL) / (sizeRH + sizeRL);
    solution.proportionRH = sum(solution.boolRxnSolution(methodOptions.RHindex)) / length(methodOptions.RHindex);
    solution.proportionRL = sum(solution.boolRxnSolution(methodOptions.RLindex) == 0) / length(methodOptions.RLindex);
    solution.proportionRxn = sum(solution.boolRxnSolution) / length(model.rxns);
end