function results = loadResultsFromZip(zipFile, pattern)
    tmp = tempname();
    mkdir(tmp);
    unzip(zipFile, tmp);
    results = {};
    files = dir([tmp filesep pattern]);
    for i = 1:numel(files)
        f = [files(i).folder filesep files(i).name];
        s = load(f, 'store');
        results{i} = s;
    end
end