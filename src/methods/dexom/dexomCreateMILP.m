function MILPproblem = dexomCreateMILP(model, options)
    % Creates a Dexom MILP problem. Code is based on the iMAT
    % implementation from COBRA Toolbox 3.0.6 with all the bells and
    % whistles for enumeration, lexicographical optimization and extra
    % parameters and constraints
    %
    % USAGE:
    %
    %    result = dexomCreateMILP(model, options)
    %
    % INPUTS:
    %    model:             input model (COBRA model structure)
    %    options:           struct with the arguments to pass to the algorithm:
    %                       * RHindex - indexes of the reactions
    %                       considered as highly expressed
    %                       * RLindex - indexes of the reactions
    %                       considered as lowly expressed
    %                       * epsilon - minimum flux for reactions
    %                       considered to be active.
    %                       * exclude - matrix. Default []
    %                       * numRL - force the solver to find a solution
    %                       with the number of lowly expressed reactions
    %                       carrying no flux indicated in the variable (default
    %                       0)
    %                       * senseRL - direction of the constraint, 'G'
    %                       indicates greater or equal to numRL, 'L' indicates
    %                       less or equal, 'E' indicates equal (default 'G')
    %                       * numRH - same as numRL but for highly
    %                       expressed reactions (default 0)
    %                       * senseRH - same as senseRL but for highly
    %                       expressed reactions (default 'G')
    %                       * maximizeDistanceToSolution - if provided,
    %                       instead of maximizing the RH+RL, it maximizes
    %                       the distance to the provided solution (as
    %                       the original [RH-fwd|RL|RH-back] binary format)
    %                       * fitScore - force an objective score
    %
    % OUTPUTS:
    %                       * MILPproblem - specification of the MILP problem to
    %                       solve to be passed to a concrete solver (use
    %                       SolveCobraMILP from Cobra for example)
    %                       * options - options used to build the MILP problem,
    %                       including the user specific options and the default
    %                       options selected by the algorithm. This is useful
    %                       also for reproducibility, as slightly different
    %                       options can produce different results.
    %
    
    if ~isfield(options,'RHindex')
        error('Indexes of the highly expressed reactions not provided\n');
    end

    if ~isfield(options,'RLindex')
        error('Indexes of the lowly expressed reactions not provided\n');
    end

    RHindex = options.RHindex;
    RLindex = options.RLindex;
    epsilon = options.epsilon;
    senseRH = options.senseRH;
    senseRL = options.senseRL;
    numRH = options.numRH;
    numRL = options.numRL;
    S = model.S;
    lb = model.lb;
    ub = model.ub;

    if isempty(RHindex)
        senseRH = [];
        numRH = [];
    end
    if isempty(RLindex)
        senseRL = [];
        numRL = [];
    end

    % Create the A matrix (constraint equations x variables) Ax <=> b
    %
    % Same A matrix as in the original iMAT code +2 rows, one for the numRL
    % constraint and other for the numRH constraint. Note that if numRL = 0,
    % numRH = 0, senseRL = 'G' and senseRH = 'G' (default values) the solution
    % of the problem is equivalent to the original iMAT.
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex), size(S,2)+2*length(RHindex)+length(RLindex));
    [m,n,s] = find(S);
    for i = 1:length(m)
        A(m(i),n(i)) = s(i);
    end

    for i = 1:length(RHindex)
        A(i+size(S,1),RHindex(i)) = 1;
        A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon;
        A(i+size(S,1)+length(RHindex),RHindex(i)) = 1;
        A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon;
    end

    for i = 1:length(RLindex)
        lbrl = lb(RLindex(i));
        ubrl = ub(RLindex(i));
        A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1;
        A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lbrl;
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1;
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ubrl;
    end

    % Constraint of the RLindex reactions to count the number of RL
    % non included in the model to force the solver to search a solution
    % with less, more or equal number of RL / RH
    lastRow = size(A, 1);
    
    for i = 1:length(RHindex)
        A(lastRow+1,i+size(S,2)+length(RHindex)+length(RLindex)) = 1;
        A(lastRow+1,i+size(S,2)) = 1;
    end

    lastRow = size(A, 1);
    for i = 1:length(RLindex)
        A(lastRow+1, i+size(S,2)+length(RHindex)) = 1;
    end
    
    % Alternatively, use fitScore to put a min constraint on the RH+RL
    % score. This can be used to change the objective function and still
    % forcing to find solutions above a threshold based on coverage of RH
    % and RL (sum all the RH/RL selected)
    A(end + 1, (size(S,2) + 1):end) = 1;

    MILPproblem.fitScoreIdx = size(A, 1);
    % Exclude all the solutions in prevSols. To do so add a new row r in
    % the A matrix per solution. For each solution s in prevSols where 
    % s == 1, add a 1 in r (r(s==1) = 1), and for each 0 in s add a -1
    % so when we perform A*x we compute the difference between current
    % selected solution in x, and previous solution s.

    % options.exclude has size N x R, where rows are vectors of presence/absence
    % of reactions in previous solutions.

    %                        2 * size(RH) + RL
    %                |-------------------------------|
    %        |-- S --|----RHf---|----RL----|---RHb---|
    %        |R1...Rn,RH1...RHm, RL1,...RLo,RH1...RHm
    %A(end,:)|0,...,0,[       constraint_vector      ]

    prevSolBounds = [];
    for i = 1:size(options.exclude, 1)
        cvector = round(options.exclude(i, :));
        % Change the RH-f/RL/RH-b binary value from the solution so each 0 is
        % -1 in order to perform the multiplication with X to count the matches
        % between the solution to exclude and the current one
        cvector(cvector == 0) = -1;
        A(end + 1, (size(S, 2)+1):end) = cvector;
        prevSolBounds(i) = sum(cvector == 1) - 1;
    end

    

    % Creating csense
    csense1(1:size(S,1)) = 'E';
    csense2(1:length(RHindex)) = 'G';
    csense3(1:length(RHindex)) = 'L';
    csense4(1:length(RLindex)) = 'G';
    csense5(1:length(RLindex)) = 'L';
    csenseFitScore = 'G';
    csense6(1:size(options.exclude,1)) = 'L'; % All less or eq for options.exclude
    csense = [csense1 csense2 csense3 csense4 csense5 senseRH senseRL csenseFitScore csense6];
    

    % Creating lb and ub
    lb_y = zeros(2*length(RHindex)+length(RLindex),1);
    ub_y = ones(2*length(RHindex)+length(RLindex),1);
    lb = [lb;lb_y];
    ub = [ub;ub_y];

    % Creating c
    c_v = zeros(size(S,2),1);
    % For the binary variables, put ones in all to count them in the
    % objective function (if the goal is to maximize coverage of RH+RL)
    c_y = ones(2*length(RHindex)+length(RLindex),1);
    
    % Add custom weights
    if ~isempty(options.rhWeights)
        c_y(1:length(RHindex)) = options.rhWeights;
        c_y((length(RHindex) + length(RLindex) + 1):end) = options.rhWeights;
    end
    if ~isempty(options.rlWeights)
        c_y(length(RHindex)+1:length(RHindex)+length(RLindex)) = options.rlWeights;
    end
    
    % Maximize differences between the reference solution and the candidate
    % by calculating the hamming distance (since the obj is c'x, we need to
    % set c based on the reference solution so c'x returns the total
    % matches between the refsol and the candidate). 
    if isfield(options, 'secondaryObjectiveWeights') && ...
        ~isempty(options.secondaryObjectiveWeights) && ...
        options.useSecondaryObjective == 1
    
        if length(options.secondaryObjectiveWeights) ~= length(c_y)
            error('Invalid referenceSolution in method options');
        end
        c_y = options.secondaryObjectiveWeights';
        % Minimize overlapping
        MILPproblem.osense = options.secondaryObjectiveSense;
        if options.fitScore <= 0
            error('Using a secondary objective for optimization without fixing the value of the first objective using fitScore');
        end
        % TODO: Instead of working directly with the MILP structure, use
        % an intermediary structure with the options required to pass
        % around to check status
        MILPproblem.secondaryObjectiveActive = 1;
    else
        % Default case: maximize gene expression agreement
        MILPproblem.osense = -1;
        MILPproblem.secondaryObjectiveActive = 0;
    end
        
    % Now c'x returns the number of coincidences between the reference
    % solution and the candidate, to maximize the distance we need to
    % minimize the number of differences
    
    c = [c_v;c_y];

    % Creating b
    b_s = zeros(size(S,1),1);
    lb_rh = lb(RHindex);
    ub_rh = ub(RHindex);
    lb_rl = lb(RLindex);
    ub_rl = ub(RLindex);
    % Calculate the bound for each prevSolution, which should be the number
    % of RH reactions with a 1 in the solution
    b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl;numRH;numRL;options.fitScore;prevSolBounds'];


    % Creating vartype
    vartype1(1:size(S,2),1) = 'C';
    vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
    vartype = [vartype1;vartype2];

    MILPproblem.A = A;
    MILPproblem.b = b;
    MILPproblem.c = c;
    MILPproblem.lb = lb;
    MILPproblem.ub = ub;
    MILPproblem.csense = csense;
    MILPproblem.vartype = vartype;
    MILPproblem.x0 = [];

end

