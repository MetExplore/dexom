function growthRatios = multiNetKnockout(model, networks, createUnionNetwork, verbose)
    if createUnionNetwork
        union = max(networks,[],1);
        nRxn = model.rxns(union == 0);
        m = removeRxns(model, nRxn);
        % Verify if there is any growth before KOs
        optimalObj = optimizeCbModel(m);
        if optimalObj.f <= 0
            warning('Base model before KO simulation have 0 flux through the biomass reaction');
            growthRatios = -1*ones(1,length(model.genes));
        else
            growthRatios = singleGeneDeletion(m)';
        end
    else
        for i = 1:size(networks, 1)
            % For each optimal network, create a new model and test
            net = networks(i,:);
            nRxn = model.rxns(net == 0);
            m = removeRxns(model, nRxn);
            % Verify if there is any growth before KOs
            optimalObj = optimizeCbModel(m);
            if optimalObj.f <= 0
                warning('Base model before KO simulation have 0 flux through the biomass reaction');
                growthRatios(i,:) = -1*ones(1,length(model.genes));
            else
                if verbose == 1
                    fprintf('Running knockouts for network %d/%d (%d reactions)...\n', i, size(networks,1), length(m.rxns));
                end
                growthRatios(i,:) = singleGeneDeletion(m);    
            end
            
        end
    end
end