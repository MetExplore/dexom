function MILPproblem = addIntegerCut(model, MILPproblem, yRef)
% Adds an integer-cut to exclude the vector yref from the original
% MILPproblem. The integer-cut is a linear constraint that prevents
% obtaining the yref vector as a solution, excluding it from the set of
% feasible solutions.
    if ~isempty(yRef) && sum(yRef) > 0 
        % Make sure they are ones and zeros, avoid precision issues
        cvector = round(yRef);
        % Change the RH-f/RL/RH-b binary value from the solution so each 0 is
        % -1 in order to perform the multiplication with X to count the matches
        % between the solution to exclude and the current one
        cvector(cvector == 0) = -1;
        % Add at the end of the list of constraints
        MILPproblem.A(end + 1, (size(model.S, 2)+1):end) = cvector;
        newCBound = sum(cvector == 1) - 1;
        MILPproblem.csense = [MILPproblem.csense 'L'];
        MILPproblem.b = [MILPproblem.b; newCBound];
    end
end