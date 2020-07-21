function evaluation = evaluateMetaEnsemble(model, zipFile, pattern, koDatasetFile)
    c = loadEssentialGenePredFromFile(zipFile, pattern);
    ko = readtable(koDatasetFile);
    ck = c == 0;
    evaluation = evaluateKO(model, ko, ck', 0.01);
end