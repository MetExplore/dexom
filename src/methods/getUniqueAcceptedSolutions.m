function us = getUniqueAcceptedSolutions(result, useIntVarFromMILP)
    if ~exist('useIntVarFromMILP','var')
        useIntVarFromMILP = 0;
    end
    if useIntVarFromMILP
        us = unique(result.intSolutions(result.intAccepted == 1,:), 'rows');
    else
        us = unique(result.solutions(result.accepted == 1,:), 'rows');
    end
end