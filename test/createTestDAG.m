function test = createTestDAG(layers, metabolitesPerLayer, method, verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    model = dagNet(layers, metabolitesPerLayer);
    numOptimalSolutions = metabolitesPerLayer^layers;
    e.maxIterations = numOptimalSolutions * 3;
    e.maxSolutions = round(numOptimalSolutions * 1.2);
    e.maxUniqueSolutions = round(numOptimalSolutions * 1.2);
    e.optTol = 1e-8;
    e.maxEnumTime = 60;
    m.RLindex = 1:length(model.rxns);
    m.RHindex = [];
    e.verbose = verbose;
    test.options = setupMethodOptions(method, m, e);
    test.model = model;
end