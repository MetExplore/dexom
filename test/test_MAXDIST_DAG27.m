function pass = test_MAXDIST_DAG27(verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    test = createTestDAG(3, 3, 'dexom-maxdist', verbose);
    result = sequentialNetworkEnumeration(test.model, test.options);
    us = getUniqueAcceptedSolutions(result);
    passSize = size(us, 1) == 27;
    % An optimal solution has 6 selected reactions and 20 non selected
    % reactions (thus an score of 20 as all reactions are in the RL set).
    % There are 2 reactions (in, out) shared by all the solutions. Two
    % maximal distant solutions cannot share more than 2 reactions, 
    % or 20-2 = 18 non selected reactions.
    x = result.intSolutions(1,:);
    y = result.intSolutions(2,:);
    passLength = sum(1-x) == 6 && sum(1-y) == 6;
    passDiff = sum(x == y) == 18;
    pass = passSize && passLength && passDiff;
end