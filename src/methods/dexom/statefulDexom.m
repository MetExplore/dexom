function result = statefulDexom(model, ctx)
% Iterative version of DEXOM. Do not run this method directly. This method
% is used by sequentialNetworkEnumeration.m. In order to run DEXOM, 
% launch the method with dexom.m.

    methodOptions = ctx.options.method;
    % Update the behavior during the main loop from the method
    % sequentialNetworkEnumeration. Variable ctx contains the context with
    % the information generated outside.
    % TODO: Check if Matlab OOP is worth it to improve the code
    if isempty(fieldnames(ctx.lastSol))
        % Initialize a new MILP using the provided settings
        milp = dexomCreateMILP(model, methodOptions);
    else
        milp = ctx.lastSol.internal.problemDef;
    end
    idxAccepted = find(ctx.accepted == 1);
    ub = milp.ub;
    lb = milp.lb;
    
    % Check if there is any valid proposition for rxn-cut. If so, modify
    % the milp problem
    if existProposedRxnCut(ctx)
        % Change bounds and restore at the end
        milp = updateReactionBounds(milp, ctx.rxnIds, ctx.constraints, methodOptions.tol); 
    end
        
    % If a secondary objective is enabled and there are no propositions
    % for reaction cuts, use them
    if methodOptions.useSecondaryObjective == 1 && ...
            ~existProposedRxnCut(ctx) && ~isempty(fieldnames(ctx.lastSol))
        % Check if the MILP problem is already configured. If it's not,
        % build a new MILP problem setting the secondary objective to
        % minimize overlap with the reference solution
        milp = checkAndSetupSecondObjective(milp, model, ctx);
        % Get last optimal solution
        if length(idxAccepted) > 1
            i = idxAccepted(end);
        else
            i = idxAccepted;
        end
        % Apply the corresponding strategy
        if strcmp(methodOptions.secondaryObjectiveMode, 'maxdist')
            c_y = ctx.intSolutions(idxAccepted(end),:)';
        elseif strcmp(methodOptions.secondaryObjectiveMode, 'dexom')
            % Update the reference solution and the overlapping set based
            % on the strategy
            t = methodOptions.distAnnealing ^ max(0, (ctx.iter - ctx.lastIterRxnCutProposal));
            % Apply only on the positions that have a one
            sel = find(ctx.intSolutions(i,:) == 1);
            rc = rand(1, length(sel)) > t;
            % Selected positions will take ones, 0s for the rest
            vec = zeros(1, size(ctx.intSolutions,2));
            vec(sel(rc)) = 1;
            c_y = vec';  
        elseif strcmp(methodOptions.secondaryObjectiveMode, 'random')
            % Get the ones of the last optimal solution
            sel = find(ctx.intSolutions(i,:) == 1);
            % Decide a size N at random and pick N ones
            n = randi([1 length(sel)]);
            rndIdx = randperm(length(sel));
            vec = zeros(1, size(ctx.intSolutions,2));
            vec(sel(rndIdx(1:n))) = 1;
            c_y = vec';  
        end
        % Redefine the optimization function
        c_v = zeros(size(model.S,2),1);
        milp.c = [c_v;c_y];    
    end
    
    % Solve the MILP problem
    solution = dexomSolveCobraMILP(milp, methodOptions);
    result.fluxSolution = solution.cont';
    if length(result.fluxSolution) ~= length(model.rxns)
        % Invalid solution returned by the solver
        result.fluxSolution = zeros(1, length(model.rxns));
        result.isValid = 0;
        if methodOptions.decreaseFitScoreAfterTries > 0 && ...
                ctx.consecutiveRejections >= methodOptions.decreaseFitScoreAfterTries
            % Get the current used fitScore
            currFitScore = milp.b(milp.fitScoreIdx);
            if currFitScore > (ctx.defSol.internal.solverOutput.obj * ...
                    (1 - ctx.options.enum.optTol)) - 1
                milp.b(milp.fitScoreIdx) = currFitScore - 1;
            end
        end
    else
        result.isValid = 1;
        % If integer cuts are enabled, exclude last solution
        if methodOptions.excludeSolutionsByIntegerCuts == 1
            % Add new integer cut using the RH/RL binary variables
            milp = addIntegerCut(model, milp, solution.int'); 
        end
    end   
    result.boolRxnSolution = abs(result.fluxSolution) >= methodOptions.tol;
    
    % Reset MILP bounds 
    milp.ub = ub;
    milp.lb = lb;
        
    result.internal.solverOutput = solution;
    result.internal.problemDef = milp;
    result.internal.problemType = 'MILP';
end

function isValid = existProposedRxnCut(ctx)
    invalid = isempty(ctx.rxnIds) || sum(ctx.rxnIds <= 0) > 0;
    isValid = ~invalid;
end

function milp = checkAndSetupSecondObjective(milp, model, ctx)
    % Check if the MILP is configured to use a second objective
    % TODO: Do not rely on MILP osense to check if the objective
    % is the right one as is not explicit enough and can lead to
    % problems in the future
    methodOptions = ctx.options.method;
    idxAccepted = find(ctx.accepted == 1);
    if isempty(methodOptions.secondaryObjectiveWeights) && ...
            milp.secondaryObjectiveActive == 0
        yRef = ctx.intSolutions(idxAccepted(end),:);
        opts = methodOptions;
        opts.secondaryObjectiveWeights = yRef;
        % If the user did not provide any custom score, use the optimal.
        if methodOptions.fitScore <= 0
            opts.fitScore = sum(opts.secondaryObjectiveWeights);
        end
        % Create a new MILP using the secondary objective initialized 
        % with the last solution accepted. This changes also the
        % optimization goal (from maximization to minimization)
        us = unique(ctx.intSolutions(ctx.intAccepted == 1,:), 'rows');
        opts.exclude = us;
        milp = dexomCreateMILP(model, opts);
    end
end