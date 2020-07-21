function MILPproblem = changeRxnBoundsMILP(model, MILPproblem, RHindex, ...
    RLindex, rxnID, lowerBound, upperBound, epsilon)
    % Directly manipulate the MILP problem th change bounds without
    % reconstructing the full matrix
    rhi = find(RHindex == rxnID);
    rli = find(RLindex == rxnID);
    if isempty(rhi), rhi = 0; end
    if isempty(rli), rli = 0; end
        
    if rhi && rli
        error('Reaction cannot be in RHindex and RLindex at the same time');
    end
    
    if rhi
        i = rhi;
        MILPproblem.A(i+size(model.S,1),rxnID) = 1;
        MILPproblem.A(i+size(model.S,1),i+size(model.S,2)) = lowerBound - epsilon;
        MILPproblem.A(i+size(model.S,1)+length(RHindex),rxnID) = 1;
        MILPproblem.A(i+size(model.S,1)+length(RHindex),i+size(model.S,2)+length(RHindex)+length(RLindex)) = upperBound + epsilon;
    end

    if rli
        i = rli;
        MILPproblem.A(i+size(model.S,1)+2*length(RHindex),rxnID) = 1;
        MILPproblem.A(i+size(model.S,1)+2*length(RHindex),i+size(model.S,2)+length(RHindex)) = lowerBound;
        MILPproblem.A(i+size(model.S,1)+2*length(RHindex)+length(RLindex),rxnID) = 1;
        MILPproblem.A(i+size(model.S,1)+2*length(RHindex)+length(RLindex),i+size(model.S,2)+length(RHindex)) = upperBound;
    end
    
    MILPproblem.lb(rxnID) = lowerBound;
    MILPproblem.ub(rxnID) = upperBound;
    
    % Find the position of the bound for the reaction. Note that b is
    % defined as [b_s;lb_rh;ub_rh;lb_rl;ub_rl;numRH;numRL;prevSolBounds'];
    % b_s is size(S, 1)
    offset = size(model.S,1);
    % lb_rh is size = length(RHindex)
    % ub_rh is size = length(RHindex)
    % lb_rl is size = length(RLindex)
    % ub_rl is size = length(RLindex)
    if rhi
        % Change bounds in lb_rh and ub_rh
        MILPproblem.b(offset + rhi) = lowerBound;
        MILPproblem.b(offset + length(RHindex) + rhi) = upperBound;
    end
    if rli
        MILPproblem.b(offset + 2*length(RHindex) + rli) = lowerBound;
        MILPproblem.b(offset + 2*length(RHindex) + length(RLindex) + rli) = upperBound;
    end
end