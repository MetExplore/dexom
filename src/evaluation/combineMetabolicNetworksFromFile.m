function [networks, methods] = combineMetabolicNetworksFromFile(zipFile, pattern)
    r = loadResultsFromZip(zipFile, pattern);
    s = [];
    methods = [];
    for i = 1:numel(r)
        result = r{i}.store.result;
        method = repelem(string(r{i}.store.info.method), sum(result.accepted == 1));
        methods = [methods, method];
        s = [s; result.solutions(result.accepted == 1, :)];
        s = logical(round(s));
    end
    networks = s;
    %[networks, ia, ic] = unique(s, 'rows');
    %methods = methods(ia);
end