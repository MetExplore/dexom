function dexomInit(cobraToolboxInitMode, verbose)
% Initializes DEXOM for diversity-based enumeration of optimal
% context-specific metabolic networks. It requires COBRA Toolbox.
% If COBRA is not in the path, it initializes the embedded version
% of COBRA Toolbox (this is the recommended setting).

    global SOLVERS;
    global CBT_MILP_SOLVER;
    global ENV_VARS;
    
    printLevel = ENV_VARS.printLevel;
    
    if ~exist('cobraToolboxInitMode','var') || isempty(cobraToolboxInitMode)
        cobraToolboxInitMode = 0;
    end
    if ~exist('verbose','var')
       verbose = 0; 
    end
    ENV_VARS.printLevel = verbose;
    
    ver = dexomVersion();
    removedPath = [];
    fprintf('\nDEXOM: Diversity-based Extraction of Optimal Metabolic-networks (v%s)\n', ver);
    fprintf('This version was tested with Matlab 2015b (CPLEX v12.8), 2018a (CPLEX v12.9),\n');
    fprintf('2018b (CPLEX v12.8) and COBRA Toolbox v3.0.6 on Windows 10.\n');
    
    paperLink = 'https://doi.org/10.1101/2020.07.17.208918';
    if usejava('desktop')
        paperLink = ['<a href=\"', paperLink, '\">', paperLink, '</a>'];
    end
    fprintf(['Research paper: ', paperLink, '\n\n']);
    
    fprintf('> Initializing DEXOM library for Matlab\n');
    PROJDIR = fileparts(which('dexomInit'));
    addpath(PROJDIR);
    addpath(genpath([PROJDIR filesep 'src']));
    addpath(genpath([PROJDIR filesep 'test']));
    addpath(genpath([PROJDIR filesep 'evaluation']));
    cd(PROJDIR);
    fprintf(' + Adding external dependencies...\n');
    if cobraToolboxInitMode == -1 && isempty(SOLVERS)
        error('COBRA Toolbox not initialized in the system. Use dexomInit(0) to load the embedded COBRA Toolbox or dexomInit(1) to load an installed version of COBRA Toolbox.');
    end
    if cobraToolboxInitMode == 1
        fprintf('> Trying to initialize COBRA Toolbox ... ');
        if isempty(SOLVERS)
            try
                initCobraToolbox(false);
            catch
                fprintf('Failed (COBRA Toolbox not in the path)\n');
                cobraToolboxInitMode = 0;
            end
        else
            fprintf('Done (already loaded)\n');
        end
    end
    
    % Use the embedded version
    if cobraToolboxInitMode == 0
        % Remove any directory related to cobratoolbox in the path
        fprintf(' + Checking and replacing previous COBRA Toolbox installations...');
        removedPath = removeCobraFromPath('initCobraToolbox');
        fprintf(' %d entries removed.\n', numel(removedPath));
        if numel(removedPath) > 0
            saveArrayStrings('cobrapath.old', removedPath); 
            fprintf('\t - In order to use your previous COBRA Toolbox, you need to reinstall it with initCobraToolbox\n');
            fprintf('\t - Or you can automatically restore your previous COBRA with restoreCobraToolboxPath()\n');
        end
        cd([PROJDIR filesep 'modules' filesep 'cobratoolbox'])
        fprintf(' + Initializing the embedded COBRA Toolbox (use dexomInit(0,1) to show log)...\n');
        initializeCobraToolboxLib();
        cd(PROJDIR);
    end
    
    % Set CPLEX as the main solver if it's installed
    sOk = 0;
    if SOLVERS.ibm_cplex.installed
        [sOk, sInstalled] = changeCobraSolver('ibm_cplex', 'all', 0);
        v = getCobraSolverVersion('ibm_cplex', 0);
        fprintf('> IBM CPLEX selected as the default solver (v%s)\n', v);
    elseif SOLVERS.gurobi.installed
        [sOk, sInstalled] = changeCobraSolver('gurobi', 'all', 0);
        v = getCobraSolverVersion('gurobi', 0);
        fprintf('> GUROBI selected as the default solver (v%s)\n', v);    
    else
        warning('> CPLEX or GUROBI not detected in the system. Use changeCobraSolver to manually specify the solver to be used');
    end
    
    if sOk == 0
        warning('> The selected solver is not working properly\n');
    end
    
    % Run a small battery of tests
    fprintf('> Testing DEXOM (solver %s) ... ', CBT_MILP_SOLVER);
    if runTests() == 1
        fprintf('Done.\n> DEXOM is ready to use.\n');
    else
        fprintf('DEXOM test(s) not passed\n');
        warning('There is a problem with your setup');
    end
    if numel(removedPath) > 0
        fprintf('\nNOTE: Your previous installation of COBRA Toolbox\n');
        fprintf('has been replaced with the embedded COBRA Toolbox v3.0.6.\n');
        fprintf('If you want to replace it with your previous COBRA,\n');
        fprintf('just call the restoreCobraToolboxPath function.\n');  
    end
    ENV_VARS.printLevel = printLevel;
end