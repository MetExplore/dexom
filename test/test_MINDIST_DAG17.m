function pass = test_MINDIST_DAG17(verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    test = createTestDAG(2, 3, 'dexom-mindist', verbose);
    result = sequentialNetworkEnumeration(test.model, test.options);
    us = getUniqueAcceptedSolutions(result);
    passSize = size(us, 1) == 9;
    % Optimal solution is a path of only 5 reactions. Any change in one
    % reaction implies changing at least another one in the last layer.
    % Thus, the minimum distance is sharing 3 reactions in the path.
    x = 1 - result.intSolutions(1,:);
    y = 1 - result.intSolutions(2,:);
    passLength = sum(x) == 5 && sum(y) == 5;
    passDiff = sum(x & y) == 3;
    pass = passSize && passLength && passDiff;
end