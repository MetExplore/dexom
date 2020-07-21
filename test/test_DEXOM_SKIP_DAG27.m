function pass = test_DEXOM_SKIP_DAG27(verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    test = createTestDAG(3, 3, 'dexom-default/norand', verbose);
    test.options.enum.enumGreedySkipReactions = 0;
    result = sequentialNetworkEnumeration(test.model, test.options);
    us = getUniqueAcceptedSolutions(result);
    passSize = size(us, 1) == 27;
    passSwitch = result.lastIterRxnCutProposal == length(test.model.rxns) + 2;
    pass = passSize && passSwitch;
end