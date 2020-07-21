function classification = evaluateKO(model, koTable, growthRatios, koTol)
    listGenes = koTable.GENE_ID;
    geneIds = findGeneIDs(model, listGenes);
    % Exclude genes for which there is no match
    notInModel = find(geneIds <= 0);
    idx = setdiff(1:length(listGenes), notInModel);
    orderedGrowthRatios = growthRatios(geneIds(geneIds > 0), 1);
    % Essential genes for Aerobic condition in Rintala Medium
    trueEssential = koTable(idx,:).LABEL_AEROBIC == 1;
    trueNonEssential = ~trueEssential;
    
    % True Positive = Essential gene with a KO growth ratio prediction <= koTol
    % True Negative = Non essential gene with a KO growth ratio predicted
    % >= koTol
    predictedEssential = orderedGrowthRatios <= koTol;
    predictedNonEssential = ~predictedEssential;
    
    tp = predictedEssential & trueEssential;
    tn = predictedNonEssential & trueNonEssential;
    fp = predictedEssential & trueNonEssential;
    fn = predictedNonEssential & trueEssential;
    
    n_tp = sum(tp);
    n_tn = sum(tn);
    n_fp = sum(fp);
    n_fn = sum(fn);
        
    classification.sensitivity = (n_tp/(n_tp+n_fn));
    classification.specificity = (n_tn/(n_tn+n_fp));
    classification.positivePredictive = (n_tp/(n_tp+n_fp));
    classification.negativePredictive = (n_tn/(n_fn+n_tn));
    classification.accuracy = (n_tp + n_tn) / (n_tp + n_fp + n_tn + n_fn);
    classification.f1 = 2 * (classification.sensitivity * classification.specificity) /...
        (classification.sensitivity + classification.specificity);
    classification.mcc = (n_tp * n_tn - n_fp * n_fn)/...
        sqrt((n_tp + n_fp)*(n_tp + n_fn)*(n_tn + n_fp)*(n_tn + n_fn));
    
    essential = zeros(1, length(model.genes));
    essential(geneIds(geneIds > 0)) = predictedEssential;
    classification.essentialGeneFromModel = sparse(essential);
    classification.N_TP = n_tp;
    classification.N_TN = n_tn;
    classification.N_FP = n_fp;
    classification.N_FN = n_fn;
    
    classification.TPR = classification.sensitivity;
    classification.FPR = 1 - classification.specificity;
    
    %classification.TP = model.genes(tp);
    %classification.TN = model.genes(tn);
    %classification.FP = model.genes(fp);
    %classification.FN = model.genes(fn);
end

