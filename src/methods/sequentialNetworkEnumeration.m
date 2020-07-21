function ctx = sequentialNetworkEnumeration(model, options)
    startTime = tic();
    lastTimeUpdate = startTime;
    ctx = initializeContext(model);
    ctx.options = options;
    
    % Run the default problem
    ctx.defSol = doReconstruction(model, ctx);
    if ctx.defSol.isValid == 0
        error('Cannot find an initial solution');
    end
    ctx.elapsedTime(ctx.iter) = toc(startTime);
    ctx = registerSolution(ctx.defSol, ctx);
    ctx.medoid = ctx.defSol.internal.solverOutput.int';
    printStatus(ctx);
    
    distances = 0;
    uniqueSolutions = ctx.defSol.activity;
    
    if size(uniqueSolutions, 1) ~= 1
        error('Expected row-vector solution');
    end
    
    while 1
        % Calculate proportion of time covered and proportion of iterations
        % covered and return the maximum
        ctx.iter = ctx.iter + 1;
        
        % Compute progress by taking into account all exit conditions
        pSols = ctx.metrics.numUniqueSolutions(end) / (options.enum.maxUniqueSolutions + 1);
        pIters = ctx.iter / (options.enum.maxIterations+1);
        pTime = min(ctx.elapsedTime(end), options.enum.maxEnumTime) / options.enum.maxEnumTime;
        ctx.progress = min(1, max(pSols, max(pIters, pTime)));
        printStatus(ctx);
        
        if checkExitConditions(ctx)
            break;
        end
        
        % Modify model based on the strategy. For each reaction in rxnIds,
        % there is a corresponding constraint which indicates the
        % constraint to apply to the given reaction (1 = block, 2 = force
        % forward flux, 3 = force backward flux)
        [ctx.rxnIds, ctx.constraints] = makeRxnCutProposal(model, ctx);
        if (isempty(ctx.rxnIds) || sum(ctx.rxnIds <= 0) > 0) && ctx.lastIterRxnCutProposal == 0
            % Check where there was a switch
            ctx.lastIterRxnCutProposal = ctx.iter;
            % Terminate if there is integer-cut methods or second obj.
            % function is not defined, as no new solutions are guaranteed
            if ctx.options.method.useSecondaryObjective == 0 && ...
                ctx.options.method.excludeSolutionsByIntegerCuts == 0
                printLn('No more proposition, terminating...\n', ctx.options.enum.verbose);
				break;
            end
        end
        
        % Update the optimal activity
        ctx.status = updateRxnStatusByConstraint(ctx.status, ctx.rxnIds, ...
            ctx.constraints);
        
        % Apply the constraints to the model and generate a new model that
        % will be used as the base model for reconstruction
        if options.method.useDefaultIterativeModelConstraintUpdate == 1
            m = applyNewConstraints(model, ctx.rxnIds, ctx.constraints, options.method.tol);
        else
            % Otherwise let the method take care of the updates (this
            % allows adhoc optimizations to compute problems faster instead
            % of regenerating the entire model to modify only a few things)
            m = model;
        end
        
        % Re-run method on the new model with the new constraints
        sol = doReconstruction(m, ctx);
        ctx.elapsedTime(ctx.iter) = toc(startTime);
        ctx = registerSolution(sol, ctx); 
        
        % Count the consecutive number of invalid solutions
        if ctx.accepted(end) == 0
            ctx.consecutiveRejections = ctx.consecutiveRejections + 1;
            ctx.totalRejected = ctx.totalRejected + 1;
        else
            ctx.consecutiveRejections = 0;
        end
        
        if options.enum.metricsUpdateFrequency >= 0 && ...
                toc(lastTimeUpdate) > options.enum.metricsUpdateFrequency
            lastTimeUpdate = tic();
            % TODO: Unique solutions can be computed using the raw integer
            % solution (where vars are RH/RL reactions) or the full solution
            % using all reactions.
            uniqueSolutions = unique(ctx.solutions(ctx.accepted == 1,:), 'rows');
            %uniqueSolutions = unique(ctx.intSolutions(ctx.intAccepted == 1,:), 'rows');
            distances = calculateDistance(uniqueSolutions, options.enum.metricsDistance);
            ctx.metrics.lastUpdate = ctx.iter;
        end
        
        % Move to method update metrics?
        % Store only summary information
        % REMOVE METRICS FROM THE CONTEXT
        ctx.metrics.numUniqueSolutions(end+1) = size(uniqueSolutions, 1);
        ctx.metrics.distance.min(end+1) = min(distances);
        ctx.metrics.distance.max(end+1) = max(distances);
        ctx.metrics.distance.diff(end+1) = ...
            ctx.metrics.distance.max(end) - ctx.metrics.distance.min(end);
        ctx.metrics.distance.mean(end+1) = mean(distances);
        ctx.metrics.distance.std(end+1) = std(distances);
        ctx.metrics.distance.coeffv(end+1) = std(distances)/mean(distances);
        % Compute the average nearest neighbor distance (dispersion) and
        % also the medoid of the solution.
        Z = squareform(distances);
        totalDist = sum(Z,2);
        [~,medoidID] = min(totalDist);
        ctx.medoid = uniqueSolutions(medoidID,:);
        Z(Z==0) = Inf;
        ctx.metrics.distance.avgNearest(end+1) = mean(min(Z,[],2));
    end
    
end

function printLn(line, verbose)
    if verbose <= 0, return; end
    fprintf(line);
end

function printStatus(ctx)
    if ctx.options.enum.verbose <= 0, return; end  
    nRxn = sum(ctx.lastSol.boolRxnSolution);
    progress = ctx.progress * 100;
    propRH = ctx.lastSol.proportionRH * 100;
    propRL = ctx.lastSol.proportionRL * 100;
    propRxn = ctx.lastSol.proportionRxn * 100;
    elapsed = datestr(datenum(0,0,0,0,0,ctx.elapsedTime(end)), 'HH:MM:SS');
    intScore = sum(ctx.lastSol.internal.solverOutput.int);
    % metrics update symbol
    metricsStatus = '(=)';
    % Change for metrics
    if isfield(ctx, 'metrics')
        unique = ctx.metrics.numUniqueSolutions(end);
        meanDistance = ctx.metrics.distance.mean(end);
        %diffDistance = ctx.metrics.distance.diff(end);
        stdDistance = ctx.metrics.distance.std(end);
        avgNearest = ctx.metrics.distance.avgNearest(end);
        if ctx.iter == ctx.metrics.lastUpdate + 1
            metricsStatus = '(*)';
        end
        metricsLine = sprintf(', metrics %s: unique = %04d, %s dist = mean %08.4f / std %08.4f / avg.near %08.4f', ...
            metricsStatus, unique, ctx.options.enum.metricsDistance, meanDistance, stdDistance, avgNearest);
    else
        metricsLine = '';
    end

    fprintf('%s %s method:%s %06.2f%% %04dit [solutions = %04d%s] Rxn-cut proposal %04d (c-%d), solution: num Rxn = %04d (%06.2f%% rxn), Obj = %04d, iScore = %04d, nScore = %.4f (%06.2f%% RH, %06.2f%% RL), accepted = %d\n', ...
        datestr(datetime('now')), elapsed, ctx.options.method.name, ...
        progress, ctx.iter, sum(ctx.accepted), metricsLine, ...
        ctx.rxnIds, ctx.constraints, nRxn, propRxn, ...
        ctx.lastSol.internal.solverOutput.obj, intScore, ...
        ctx.lastSol.normalizedScore, propRH, propRL, ctx.accepted(end));
end

function ctx = registerSolution(sol, ctx)
    ctx.lastSol = sol;
    ctx.solutions(ctx.iter,:) = sol.activity;
    if sol.isValid
        ctx.intSolutions(ctx.iter,:) = sol.internal.solverOutput.int';
        ctx.contSolutions(ctx.iter,:) = sol.internal.solverOutput.cont';
    else
        ctx.intSolutions(ctx.iter,:) = zeros(1, 2 * length(ctx.options.method.RHindex) + length(ctx.options.method.RLindex));
        ctx.contSolutions(ctx.iter,:) = zeros(1, size(ctx.skipRxn, 2));
    end
    obj = sol.internal.solverOutput.obj;
    if isempty(obj)
        obj = nan;
    end    
    ctx.objectives(ctx.iter) = obj;
    ctx.scores(ctx.iter,:) = sol.normalizedScore;
    ctx.isValid(ctx.iter) = sol.isValid;
    accept = sol.isValid && sol.normalizedScore >= ...
        ctx.defSol.normalizedScore * (1 - ctx.options.enum.optTol);
    ctx.intAccepted(ctx.iter) = accept;
    % TODO: unify the intAccepted and accepted
    ctx.accepted(ctx.iter) = accept;
    if accept
        if ctx.options.enum.enumGreedySkipReactions == 1
            ctx.skipRxn = updateRxnStatus(sol.activity, ctx.skipRxn);
        end    
    end
end

function status = updateRxnStatus(activityRxnSolution, status)
    if isempty(status)
        status = zeros(3, length(activityRxnSolution));
    end
    status(1, activityRxnSolution ==  0) = 1; % Rxn not present in optimal sol
    status(2, activityRxnSolution ==  1) = 1; % Rxn present with fwd flux in optimal sol
    status(3, activityRxnSolution == -1) = 1; % Rxn present with reversed flux in optimal sol
end

function status = updateRxnStatusByConstraint(status, rxnIds, constraints)
    for i = 1:length(rxnIds)
        id = rxnIds(i);
        constraint = constraints(i);
        if id > 0 && constraint > 0
            status(constraint, id) = 1;
        end
    end
end

function [rxnId, test] = makeRxnCutProposal(model, ctx)
    switch ctx.options.enum.enumStrategy
        case 'default'
            [rxnId, test] = defaultSelection(model, ctx, 1, 1, 1);
        case 'block'
            [rxnId, test] = defaultSelection(model, ctx, 1, 0, 0);
        case 'forward'
            [rxnId, test] = defaultSelection(model, ctx, 0, 1, 0);
        case 'backward'
            [rxnId, test] = defaultSelection(model, ctx, 0, 0, 1);
        case 'random'
            [rxnId, test] = randomSelection(model, ctx, 1, 1, 1);
        otherwise
            rxnId = 0;
            test = 0;
    end
end

function available = getCandidateReactions(model, ctx)
    available(1,:) = ctx.status(1,:) == 0 & ctx.skipRxn(1,:) == 0 & ...
        model.ub' > ctx.options.method.tol & ...
        model.lb' < -ctx.options.method.tol;
    available(2,:) = ctx.status(2,:) == 0 & ctx.skipRxn(2,:) == 0 &...
        model.ub' > ctx.options.method.tol;
    available(3,:) = ctx.status(3,:) == 0 & ctx.skipRxn(3,:) == 0 & ...
        model.lb' < -ctx.options.method.tol;
end

function [rxnId, test] = defaultSelection(model, ctx, block, forward, backward)
    available = getCandidateReactions(model, ctx);
    if block
        blck = min(find(available(1,:)==1));
        if isempty(blck)
            blck = inf;
        end
    else
        blck = inf;
    end
    if forward
        fwd = min(find(available(2,:)==1));
        if isempty(fwd)
            fwd = inf;
        end
    else
        fwd = inf;
    end
    if backward
        bckw = min(find(available(3,:)==1));
        if isempty(bckw)
            bckw = inf;
        end
    else
        bckw = inf;
    end
        
    rxnId = min(blck, min(fwd, bckw));
    if rxnId == blck
        test = 1;
    elseif rxnId == fwd
        test = 2;
    elseif backward
        test = 3;
    else
        error('Invalid option');
    end 
    
    if isempty(rxnId) || rxnId > length(model.rxns)
        rxnId = 0;
        test = 0;
    end
end


function [rxnId, test] = randomSelection(model, ctx, block, forward, backward)
    available = getCandidateReactions(model, ctx);        
    test = [];
    rxnId = 0;
    
    if sum(available(1,:)) > 1 && block, test = [test 1]; end
    if sum(available(2,:)) > 1 && forward, test = [test 2]; end
    if sum(available(3,:)) > 1 && backward, test = [test 3]; end
    
    if length(test) > 1
        test = randsample(test, 1);
    end
    
    rxnId = find(available(test,:));
    if length(rxnId) > 1
        rxnId = randsample(rxnId, 1);
    end
    
end

function m = applyNewConstraints(model, rxnIds, constraints, tol)
    m = model;
    for i = 1:length(rxnIds)
        id = rxnIds(i);
        constraint = constraints(i);
        if id > 0 && constraint > 0
            switch constraint
                case 1 % Block reaction
                    m.lb(id) = 0;
                    m.ub(id) = 0;
                case 2 % Force forward flux
                    m.lb(id) = 2*tol;
                case 3 % Force backward flux
                    m.ub(id) = -2*tol;
            end
        end
    end
end

function distances = calculateDistance(solutions, method)
    if size(solutions, 2) > 1
        distances = pdist(solutions, method);
        if isempty(distances) || sum(isnan(distances)) > 0
            distances = 0;
        end
    else
        distances = 0;
    end
end



function ctx = initializeContext(model)
    ctx.iter = 1;
    ctx.solutions = [];
    ctx.contSolutions = [];
    ctx.intSolutions = [];
    ctx.scores = [];
    ctx.objectives = [];
    ctx.accepted = [];
    ctx.intAccepted = [];
    ctx.elapsedTime = [];
    ctx.progress = 0;
    ctx.uniqueLastUpdate = 0;
    ctx.rxnIds = 0;
    ctx.constraints = 0;
    ctx.isValid = [];
    ctx.status = zeros(3, length(model.rxns));
    ctx.skipRxn = zeros(3, length(model.rxns));
    ctx.lastSol = struct;
    ctx.defSol = struct;
    ctx.medoid = nan;
    ctx.consecutiveRejections = 0;
    ctx.totalRejected = 0;
    ctx.lastIterRxnCutProposal = 0;
    ctx.metrics.distance.min = 0;
    ctx.metrics.distance.max = 0;
    ctx.metrics.distance.diff = 0;
    ctx.metrics.distance.mean = 0;
    ctx.metrics.distance.std = 0;
    ctx.metrics.distance.coeffv = 0;
    ctx.metrics.distance.avgNearest = 0;
    ctx.metrics.numUniqueSolutions = 1;
    ctx.metrics.lastUpdate = 0;
end

function exit = checkExitConditions(ctx)
    exit = 1;
    v = ctx.options.enum.verbose;
    if ctx.iter >= ctx.options.enum.maxIterations
        printLn('Maximum number of iterations reached\n', v);
    elseif ctx.elapsedTime(end) > ctx.options.enum.maxEnumTime
        printLn('Maximum time reached\n', v);
    elseif ctx.options.enum.maxUniqueSolutions > 0 && ...
            ctx.metrics.numUniqueSolutions(end) >= ctx.options.enum.maxUniqueSolutions
        printLn('Reached limit on max unique solutions\n', v);
    elseif sum(ctx.intAccepted) > ctx.options.enum.maxSolutions
        printLn('Reached limit on max accepted solutions\n', v);
    elseif ctx.consecutiveRejections >= ctx.options.enum.maxTries && ...
            (isempty(ctx.rxnIds) || ctx.rxnIds == 0)
        printLn('Reached limit on max retries without no new solution\n', v);
    else
        exit = 0;
    end
end