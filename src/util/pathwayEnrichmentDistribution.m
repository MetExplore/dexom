function [scores, pathways] = pathwayEnrichmentDistribution(model, r)
    scores = [];
    nSols = size(r.unique, 1);
    for i = 1:nSols
        fprintf('%d/%d\n',i,nSols);
        [~, enrichment] = pathwayEnrichment(model, model.rxns(r.unique(i,:)==1));
        pval = cell2mat(enrichment(:,end));
        score = -log10(pval);
        scores = [scores score];
    end

    figure;
    hold on;
    boxplot(scores', 'labels', enrichment(:,1));
    xtickangle(90);
    yline(-log10(0.05));
    hold off;
    pathways = enrichment(:,1);
end