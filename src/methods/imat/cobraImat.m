function result = cobraImat(model, methodOptions)
% iMAT method as implemented COBRA Toolbox v3.0.6
    RHindex = methodOptions.RHindex;
    RLindex = methodOptions.RLindex;
    core = {};
    logfile = 0;
    epsilon = methodOptions.epsilon;
    runtime = methodOptions.runtime;
    printLevel = methodOptions.printLevel;

    if isfield(model,'C') || isfield(model,'E')
        issueConfirmationWarning('iMat does not handle the additional constraints and variables defined in the model structure (fields .C and .E.)\n It will only use the stoichiometry provided.');
    end
    
    %Manually add defined core reactions to the core
    if ~isempty(core)
        for i = 1:length(core)
            rloc = find(ismember(model.rxns, core{i}));
            if ~isempty(rloc) && isempty(intersect(RHindex,rloc))
                RHindex(end+1) = rloc;
            end
            if isempty(rloc)
                disp(['Manual added core reaction: ', core{i}, ' not found'])
            end
        end
    end

    S = model.S;
    lb = model.lb;
    ub = model.ub;

    % Creating A matrix
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
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
        A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1;
        A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i));
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1;
        A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i));
    end

    % Creating csense
    csense1(1:size(S,1)) = 'E';
    csense2(1:length(RHindex)) = 'G';
    csense3(1:length(RHindex)) = 'L';
    csense4(1:length(RLindex)) = 'G';
    csense5(1:length(RLindex)) = 'L';
    csense = [csense1 csense2 csense3 csense4 csense5];

    % Creating lb and ub
    lb_y = zeros(2*length(RHindex)+length(RLindex),1);
    ub_y = ones(2*length(RHindex)+length(RLindex),1);
    lb = [lb;lb_y];
    ub = [ub;ub_y];

    % Creating c
    c_v = zeros(size(S,2),1);
    c_y = ones(2*length(RHindex)+length(RLindex),1);
    c = [c_v;c_y];

    % Creating b
    b_s = zeros(size(S,1),1);
    lb_rh = lb(RHindex);
    ub_rh = ub(RHindex);
    lb_rl = lb(RLindex);
    ub_rl = ub(RLindex);
    b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];

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
    MILPproblem.osense = -1;
    MILPproblem.x0 = [];

    solution = solveCobraMILP(MILPproblem, 'timeLimit', runtime, 'logFile', logfile, 'printLevel', printLevel);
    
    % Add code to return additional information useful for different steps
    result.fluxSolution = solution.cont';
    if length(result.fluxSolution) ~= length(model.rxns)
        result.fluxSolution = zeros(1, length(model.rxns));
        result.isValid = 0;
    else
        result.isValid = 1;
    end   
    
    result.internal.problemDef = MILPproblem;
    result.internal.problemType = 'MILP';
    result.internal.solverOutput = solution;
    
    
    result.boolRxnSolution = abs(result.fluxSolution) >= methodOptions.tol;
end

