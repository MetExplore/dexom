function model = small4S(printLevel)
    if ~exist('printLevel','var')
        printLevel=0;
    end
    
     ReactionFormulas = {
        'A -> B',... %R1
        'B -> C',... %R2
        'C -> D',... %R3
        'D -> H',... %R4
        'A -> E',... %R5
        'E -> F',... %R6
        'F -> G',... %R7
        'G -> H',... %R8
        'A -> H',... %R9
        '-> 2 A',... %R10
        'H ->'       %R11
    };

    ReactionNames = {
        'R1', 'R2', 'R3', 'R4', 'R5', 'R6',...
        'R7', 'R8', 'R9', 'R10', 'R11'
    };
 
    lb = zeros(1, 11);
    ub = ones(1, 11)*1000;

    model = createModel(ReactionNames, ReactionNames,...
        ReactionFormulas, 'lowerBoundList', lb, 'upperBoundList', ub,...
        'printLevel', printLevel);
    
    model = changeObjective(model, 'R11');

    % Create the expression data for iMAT
    model.options.epsilon = 0.01;
    model.options.tol = 1e-6;
    
    % Reaction exprs in the format -1=low, 0=unknown, 1=high
    % Example to generate a set of optimal solutions with different
    % number of RH and RL
    
    model.options.reaction_levels = ...
        [0; 1; -1; 0; 0; 1; -1; 0; 1; 0; 0];
        %R1 R2 R3  R4 R5 R6 R7 R8 R9 R10 R11
        
    model.options.RHindex = find(model.options.reaction_levels == 1);
    model.options.RLindex = find(model.options.reaction_levels == -1);
end