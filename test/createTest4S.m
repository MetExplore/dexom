function test = createTest4S(method, verbose)
    if ~exist('verbose','var')
        verbose=0;
    end
    model = small4S();
    e.maxIterations = 20;
    e.maxSolutions = 20;
    e.maxUniqueSolutions = 20;
    e.optTol = 1e-8;
    e.maxEnumTime = 20;
    m.RLindex = model.options.RLindex;
    m.RHindex = model.options.RHindex;
    e.verbose = verbose;
    test.options = setupMethodOptions(method, m, e);
    test.model = model;
end