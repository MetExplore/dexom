function solution = dexomSolveCobraMILP(MILPproblem, options)
    global CBT_MILP_SOLVER;
    % Add solver-specific configuration
    if options.useRandomSeed, options.rseed = randi(intmax('int16')); end
    if ~isfield(options,'solver')
        options.solver = struct;
    end
    if strcmp(CBT_MILP_SOLVER, 'ibm_cplex')
        options.solver.output.clonelog = -1;
        options.solver.workdir = tempdir();
        if ~isnan(options.rseed)
            options.solver.randomseed = options.rseed;
        end
    elseif strcmp(CBT_MILP_SOLVER, 'gurobi')
        if ~isnan(options.rseed)
            options.solver.Seed = options.rseed;
        end
    elseif isempty(CBT_MILP_SOLVER)
        error('Undefined solver');
    end
    solution = solveCobraMILP(MILPproblem, options.solver);
    % Round to avoid having precision errors in the output.
    solution.int = round(solution.int);
end