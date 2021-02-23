function model = small4M()
    ReactionFormulas = {        % Low < 3, High > 7
        '-> A',...              % EX_A (?) -1
        '-> D',...              % EX_D (?) -1
        '-> X',...              % EX_X (?) -1
        '-> C',...              % EX_C (?) -1
        '-> Z',...              % EX_Z (?) -1
        'G ->',...              % EX_G (?)  5
        'Y ->',...              % EX_Y (?) -1
        'T ->',...              % EX_T (?) -1
        'A -> B',...            % RAB  (-)  1
        'B <=> C',...           % RBC  (?) -1
        'F -> G',...            % RFG  (+)  8
        'C <=> F',...           % RCF  (-)  2
        'D -> E',...            % RDE  (-)  1___
        'E + X -> F + Y',...    % REF1 (?) -1   |
        'E -> F',...            % REF2 (?) -1   |--> redundancy of paths
        'E + Z -> F + T',...    % REF3 (?) -1 __|
        'B <=> E'               % RBE  (?) -1
    };

    ReactionNames = {
        'EX_A', 'EX_D', 'EX_X', 'EX_C', 'EX_Z', 'EX_G', 'EX_Y', 'EXT_T',...
        'RAB', 'RBC', 'RFG', 'RCF', 'RDE', 'REF1', 'REF2', 'REF3', 'RBE'
    };

    lowerbounds = zeros(1, length(ReactionFormulas));
    upperbounds = zeros(1, length(ReactionFormulas)) + 1000;

    model = createModel(ReactionNames, ReactionNames, ReactionFormulas,...
    'lowerBoundList', lowerbounds, 'upperBoundList', upperbounds);

    % Make RCF and RBE reversible
    model = changeRxnBounds(model, 'RCF', -1000, 'l');
    model = changeRxnBounds(model, 'RBE', -1000, 'l');

    % Making RBC reversible allows to obtain the unique optimal solution
    % with an adequacy of 4.
    % model = changeRxnBounds(model, 'RBC', -1000, 'l');

    % Create the expression data for iMAT
    model.imatopt.epsilon = 1;
    model.imatopt.tol = 0.01;
    % Reaction exprs in the format -1=low, 0=unknown, 1=high
    model.imatopt.rxnValues = [0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 1; -1; -1; 0; 0; 0; 0];    
end





