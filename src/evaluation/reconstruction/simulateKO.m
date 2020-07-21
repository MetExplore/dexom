function growthRatio = simulateKO(model, net)
    % For each optimal network, create a new model and test
    nRxn = model.rxns(net == 0);
    m = removeRxns(model, nRxn);
    % Verify if there is any growth before KOs
    optimalObj = optimizeCbModel(m);
    if optimalObj.f <= 0
        warning('Base model before KO simulation have 0 flux through the biomass reaction');
        growthRatio = ones(length(model.genes),1);
    else
        growthRatio = geneKO(m);
    end
    if size(growthRatio,1) ~= length(model.genes)
        if size(growthRatio,2) == length(model.genes)
            growthRatio = growthRatio';
        else
            warning('Incorrect KO vs WT growth ratio for the network');
            growthRatio = ones(length(model.genes), 1);
        end
    end
end