function counts = loadEssentialGenePredFromFile(zipFile, pattern)
    r = loadResultsFromZip(zipFile, pattern);
    counts = zeros(1, 900);
    for i = 1:numel(r)
        counts = counts + full(r{i}.store.result.evalOrKO.essentialGeneFromModel);
    end
end