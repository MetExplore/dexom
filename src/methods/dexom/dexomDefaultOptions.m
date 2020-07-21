function options = dexomDefaultOptions()
    % Epsilon is the iMAT eps parameter. It corresponds to the absolute
    % minimal flux a highly expressed reaction has to carry to be 
    % considered active. The main idea is that highly expressed reactions
    % should carry at least a representative flux to be considered "highly
    % expressed". However, once a solution has been found, other reactions
    % from the highly expressed set of reaction can carry some flux below
    % epsilon and above tolerance, but they don't score 
    % (i.e., they are present in the model (or active)
    % but they don't score for the objective function).
    options.epsilon = 1e-4;
    % Tolerance is the minimum value to consider a reaction to carry some
    % flux in order to include it in the final model. This parameter is
    % closely related to the solver's LP tolerance. There is no need to tune
    % this parameter in general.
    options.tol = 1e-6;
    % numRL and sensRL can be used to force the solver to find a solution
    % with a number of lowly expressed reactions that is greater or equal
    % numRL (senseRL='G'), equal to numRL (senseRL='E') or less or equal
    % than numRL (senseRL='L')
    options.numRL=0;
    options.senseRL='G';
    % Same as numRL/senseRL but for the highly expressed reactions. Note
    % that forcing a concrete solution increases the complexity of the
    % problem. If a solution with a concrete number of lowly expressed
    % reactions / highly expressed reactions does not exists, it can take
    % lot of time for the solver to get an answer.
    options.numRH=0;
    options.senseRH='G';
    % Lower bound for the objective function
    options.fitScore=0;
    % By default, if no solution is found and fitScore is above optimal
    % solution * (1 - optTol) - 1, decrease one unit after N tries
    options.decreaseFitScoreAfterTries = 3;
    % Force the solver to exclude a list of solutions (for enumerating new
    % optimal solutions). Each row should be a binary vector of size
    % 2*length(RH) + length(RL) which corresponds to a solution.int value 
    % returned by the solver (see solveCobraMILP)
    options.exclude=[];
    % Random seed for the solver (only CPLEX/GUROBI).
    % NOTE: Selecting one rseed affects the performance of CPLEX 12.8, it
    % is recommended to leave rseed as NaN unless you are interested in
    % having more chances of generating different solutions with
    % exomSampling.
    options.rseed=NaN;
    % If set to 1, generate a new options.rseed in each call
    options.useRandomSeed = 0;
    options.excludeSolutionsByIntegerCuts = 1;
    % If set to 1, change the objective function to maximize the
    % differences with a reference solution provided
    options.useSecondaryObjective = 0;
    options.secondaryObjectiveMode = 'dexom';
    options.secondaryObjectiveSense = 1;
    % Used only when performing enumeration with Dexom. The parameter
    % modifies the speed at which the jumps go from close to far over time
    options.distAnnealing = 0.995;
    options.secondaryObjectiveWeights = [];
    options.name = 'dexom';
end