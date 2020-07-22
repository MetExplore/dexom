function dexomInit(cobraToolboxInitMode)
% Initializes DEXOM for diversity-based enumeration of optimal
% context-specific metabolic networks. It requires COBRA Toolbox.
% If COBRA is not in the path, it initializes the embedded version
% of COBRA Toolbox (this is the recommended setting).

    global SOLVERS;
    
    % By default don't initialize the embedded COBRA Toolbox
    % cobraToolboxInitMode:
    % - 0: don't load cobra 
    % - 1: try to see if initCobraToolbox is in the path. If not,
    % load the embedded version
    % - 2: try to load the embedded version 
    if ~exist('cobraToolboxInitMode','var') || isempty(cobraToolboxInitMode)
        cobraToolboxInitMode = 2;
    end
    
    ver = dexomVersion();
    fprintf('\nDEXOM: Diversity-based Extraction of Optimal Metabolic-networks (v%s)\n', ver);
    fprintf('This version was tested with Matlab 2015b (CPLEX v12.8), 2018a (CPLEX v12.9), 2018b (CPLEX v12.8) and COBRA Toolbox v3.0.6 on Windows 10.\n\n');
    fprintf('> Initializing DEXOM library for Matlab\n');
    PROJDIR = fileparts(which('dexomInit'));
    addpath(PROJDIR);
    addpath(genpath([PROJDIR filesep 'src']));
    addpath(genpath([PROJDIR filesep 'test']));
    addpath(genpath([PROJDIR filesep 'evaluation']));
    cd(PROJDIR);
    fprintf(' + Adding external dependencies...\n');
    if cobraToolboxInitMode == 0 && isempty(SOLVERS)
        error('COBRA Toolbox not initialized in the system. Use dexomInit(1) to try to load the installed or the embedded COBRA Toolbox.');
    end
    if cobraToolboxInitMode == 1
        fprintf('\t> Trying to initialize COBRA Toolbox ... ');
        if isempty(SOLVERS)
            try
                initCobraToolbox(false);
            catch
                fprintf('Failed (COBRA Toolbox not in the path)\n');
                cobraToolboxInitMode = 2;
            end
        else
            fprintf('Done (already loaded)\n');
        end
    end
    
    % Use the embedded version
    if cobraToolboxInitMode == 2
        % Remove any directory related to cobratoolbox in the path
        removeCobraFromPath();
        cd([PROJDIR filesep 'modules' filesep 'cobratoolbox'])
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
    fprintf('> Testing DEXOM ... ');
    if runTests() == 1
        fprintf('Done.\n> DEXOM is ready to use.\n');
    else
        fprintf('DEXOM test(s) not passed. Check if everything is well configured.\n');
    end
end