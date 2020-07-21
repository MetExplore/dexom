function pass = test_DEXOM_4S(verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    test = createTest4S('dexom-default/norand', verbose);
    test.options.enum.enumGreedySkipReactions = 0;
    result = sequentialNetworkEnumeration(test.model, test.options);
    us = getUniqueAcceptedSolutions(result);
    pass = size(us, 1) == 4;
end