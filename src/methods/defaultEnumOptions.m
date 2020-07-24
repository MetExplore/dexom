function enumOptions = defaultEnumOptions()
% Creates a structure with default values for enumeration of optimal
% context-specific metabolic networks. The fields are:
%
% * maxIterations (default 1000). Maximum number of iterations. Each
% iteration computes a new optimal network. Depending on the configuration
% of the search, these solutions might be unique (if integer-cuts are
% enabled) or not.
% * maxSolutions (default 1000). Maximum number of solutions (duplicated
% solutions are counted as well).
% * maxUniqueSolutions (default 1000). Maximum number of unique solutions.
% Unique solutions are updated using the metricsUpdateFrequency (in
% seconds). Every X seconds, metrics are updated. At this point, if unique
% solutions are above the threshold, the search is stopped.
% * maxEnumTime (default 3600). Maximum time (in seconds) for the
% enumeration of alternative optimal solutions.
% * metricsUpdateFrequency (default 0). Update frequency of metrics 
% (in seconds), such as number of unique solutions, 
% distance between solutions, etc. A value < 0 disables the calculation of
% metrics.
% * enumStrategy (default 'default'). Enumeration strategy for the
% generation of alternative solutions. The 'default' strategy corresponds
% to the 'dexom' strategy, in which an initial set of solutions is
% calculated using the 'rxn-enum' strategy and then the set of solutions
% are progressively expanded including solutions that are more distant,
% increasing the distance iteratively until distance is maximal.
% * optTol (default 1e-4). Before enumeration, an initial solutions is
% obtained using the default method options for context-specific
% reconstruction. The optimal objective score is used as a reference.
% Alternative optimal solutions that have a difference in score below 1e-4
% are still accepted as optimal solutions.
% * enumGreedySkipReactions (default 0). If = 1, the rxn-enum method skips
% those reactions that are already part of some alternative optimal network
% * maxTries (default 6). The search for optimal solutions is subject to 
% a certain degree of randomness. Certain problems may not be easy to solve 
% in the maximum time established, making the problem infeasible. maxTries 
% indicates the maximum amount of iterations that are executed after a 
% problem is detected infeasible.
% * metricsDistance (default 'hamming'). Distance function used to compute
% pairwise distance and nearest neighbor distance in the set of optimal
% solutions. Each time a new solution is added to the set of optimal
% solutios, the distances are recalculated.
enumOptions.metricsUpdateFrequency = 0;
enumOptions.maxIterations = 1000;
enumOptions.maxEnumTime = 3600;
enumOptions.maxSolutions = 1000;
enumOptions.optTol = 1e-4;
enumOptions.enumStrategy = 'default';
enumOptions.metricsDistance = 'hamming';
enumOptions.enumGreedySkipReactions = 0;
enumOptions.verbose = 1;
enumOptions.maxTries = 6;
enumOptions.maxUniqueSolutions = 1000;
end
