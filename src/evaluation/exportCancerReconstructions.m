function t = exportCancerReconstructions(cell, threshold, methods)
    s = [];
    labels = [];
    for i = 1:numel(methods)
        fname = sprintf('%s_%d_%d_%s.mat', cell, threshold(1), threshold(2), methods{i});
        fprintf('Loading %s\n', fname);
        % module evaluation has to be in the path
        try
            load(fname)
        catch E
            fprintf('Cannot locate result file %s, make sure that the sub-module evaluation is in the matlab path\n', fname);
            rethrow(E)
        end
        s = [s; uniqueDiscretized];
        labels = [labels repelem(methods(i), size(uniqueDiscretized, 1))];
    end
    t = array2table(s);
    t.label = labels';
end