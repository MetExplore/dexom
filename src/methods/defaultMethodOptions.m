function methodOptions = defaultMethodOptions()
    % By default, tolerance is set to  the default value used in COBRA
    % for LP problems. Flux below LP tolerance cannot be considered
    % active. 
    methodOptions.tol = 1e-6;
    % Epsilon is the smallest flux required in order to consider a highly
    % expressed reaction as such.
    methodOptions.epsilon = 1e-4;
    methodOptions.printLevel = 0;
    methodOptions.runtime = 300;
    methodOptions.useDefaultIterativeModelConstraintUpdate = 1;
end