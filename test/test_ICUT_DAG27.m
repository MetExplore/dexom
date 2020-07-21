function pass = test_ICUT_DAG27(verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    test = createTestDAG(3, 3, 'dexom-icut', verbose);
    result = sequentialNetworkEnumeration(test.model, test.options);
    us = getUniqueAcceptedSolutions(result);
    pass = size(us, 1) == 27;
end