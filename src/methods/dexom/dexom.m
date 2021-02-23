function [result, usedOptions] = dexom(model, methodOptions, enumOptions)
% Runs the default DEXOM algorithm for context-specific network
% reconstruction and enumeration. The default DEXOM uses the rxn-enum
% strategy to obtain an initial set of optimal solutions from which it
% expands the search looking for alternative solutions that are more
% different. 
%
% The default objective function is the same as in iMAT, in
% which the optimal sub-networks are the ones that maximize the tradeoff
% between selecting reactions associated with highly expressed enzymes and
% removing reactions associated with lowly expressed enzymes.
%
% USAGE:
%
%    result = dexom(model, methodOptions, enumOptions)
%
% INPUTS:
%    model:             input model (COBRA model structure)
%    methodOptions:     structure with the options for context-specific
%                       network reconstruction. Mandatory fields are:
%                       - RHindex: indexes of the reactions associated with
%                       highly expressed enzymes
%                       - RLindex: indexes of the reactions associated with
%                       lowly expressed enzymes.
%                       For documentation about the other available
%                       options, see defaultMethodOptions.m
%
% OPTIONAL INPUTS:
%    enumOptions:       structure with the options for optimal network
%                       enumeration. By default, the maximum number of
%                       networks to enumerate is 1,000, with an optimal
%                       tolerance of 1e-4 (enumerated networks with an
%                       objective score above optimal * (1 - 1e-4) are
%                       accepted as optimal alternative networks.
%
% OUTPUT:
%    result:            structure with the result of the enumeration.
    if ~isfield(methodOptions,'name')
        methodOptions.name = 'dexom-default';
    end
    if ~exist('enumOptions','var'), enumOptions = struct; end 
    o = setupMethodOptions(methodOptions.name, methodOptions, enumOptions);
    result = sequentialNetworkEnumeration(model, o);
    usedOptions = o;
end

