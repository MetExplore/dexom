function model = dagNet(numLayers, numMetabolitesPerLayer, printLevel)
    if ~exist('printLevel','var')
        printLevel=0;
    end
    model = createModel();
    nl = numLayers + 2;
    for i = 1:nl
        if i == 1
            model = addReaction(model,'R_in','reactionFormula','-> Met11', 'printLevel',printLevel);
            continue;
        elseif i == nl
            for j = 1:numel(layerMets)
                prevMet = layerMets{j};
                reactName = sprintf('R_%s_SINK', prevMet);
                formula = sprintf('%s -> MetSINK', prevMet);
                model = addReaction(model,reactName,'reactionFormula',formula,'printLevel',printLevel);
            end
        else
            for j = 1:numMetabolitesPerLayer
                % Connect previous metabolites with the new ones. Add 2*M
                % reactions connecting the previous M metabolites with each
                % metabolite in this layer
                metName = sprintf('Met%d%d',i,j);
                layerMets{j} = metName;
                numMets = numMetabolitesPerLayer;
                if i == 2
                    numMets = 1;
                end
                for k = 1:numMets
                    prevMet = sprintf('Met%d%d',i-1,k);
                    reactName = sprintf('R_%s_%s', prevMet, metName);
                    formula = sprintf('%s -> %s', prevMet, metName);
                    model = addReaction(model,reactName,'reactionFormula',formula,'printLevel',printLevel);    
                end
            end
        end
    end
    model = addReaction(model,'R_out','reactionFormula','MetSINK ->','printLevel',printLevel);
    model = changeObjective(model, 'R_out');
    model.lb(1) = 1;
    model.lb(end) = 1;
end