function [listRxn, fluxes] = optimalSparseFBA(model, fluxTol, optTol)
% Use DEXOM to solve the sparse FBA problem.
    if ~exist('optTol','var')
        optTol = 0.01;
    end
    if ~exist('fluxTol','var')
        fluxTol = 1e-6;
    end
    methodOptions.RLindex = 1:length(model.rxns);
    methodOptions.RHindex = [];
    methodOptions.epsilon = fluxTol;
    methodOptions.tol = fluxTol;
    enumOptions.maxIterations = 0;
    methodOptions.solver.relMipGapTol = optTol;
    [result, ~] = dexom(model, methodOptions, enumOptions);
    fluxes = result.contSolutions(1,:);
    listRxn = model.rxns(abs(fluxes) > fluxTol);
end