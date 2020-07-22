function removed = removeCobraFromPath(initScriptName, verbose)
    if ~exist('verbose','var') || isempty(verbose)
        verbose = 0;
    end
    folders = strsplit(path, ';');
    root = fileparts(which(initScriptName));
    removed = [];
    if ~isempty(root)
        for i = 1:numel(folders)
            if contains(folders{i}, root)
                removed{end + 1} = folders{i};
                rmpath(folders{i})
                if verbose == 1
                    fprintf('Removing from path %s\n', folders{i});
                end
            end
        end
    end
end

