model = dagNet(3,3);
e.maxIterations = 300;
e.maxSolutions = 300;
%e.metricsUpdateFrequency = 1e-8;
e.maxUniqueSolutions = 250;
e.optTol = 1e-8;
%e.enumStrategy = 'integer-cut';
e.maxEnumTime = 60*60*10;
m.RLindex = 1:length(model.rxns);
m.RHindex = [];
m.verbose = 0;

o = setupMethodOptions('dexom-default', m, e);
solution = sequentialNetworkEnumeration(model, o);

%solution = exom(model, m, e);
%sols = solution.intSolutions(solution.intAccepted == 1,:);
%us = unique(sols,'rows');
us = unique(solution.solutions(solution.accepted == 1,:), 'rows');
ius = unique(solution.intSolutions(solution.intAccepted == 1,:), 'rows');

[loadings,scores,latent,tsquared,explained,mu] = pca(us);
gscatter(scores(:,1), scores(:,2))

%unique(ctx.solutions(ctx.accepted == 1,:), 'rows')