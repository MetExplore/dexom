function pass = runTests(verbose)
    if ~exist('verbose','var')
        verbose = 0;
    end
    testFiles = dir(fullfile('test', 'test*.m'));
    pass = 1;
    for i = 1:size(testFiles, 1)
        [~, name, ~] = fileparts(testFiles(i).name);
        if verbose
            fprintf('Testing %s... ', name);
        end
        pass = pass & feval(name);
        if verbose
            if pass
                fprintf('PASS\n');
            else
                fprintf('FAIL\n');
            end
        end
    end
end