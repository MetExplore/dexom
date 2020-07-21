function m = updateReactionBounds(m, rxnIds, constraints, tol)
    for i = 1:length(rxnIds)
        id = rxnIds(i);
        constraint = constraints(i);
        if id > 0 && constraint > 0
            switch constraint
                case 1 % Block reaction
                    m.lb(id) = 0;
                    m.ub(id) = 0;
                case 2 % Force forward flux
                    m.lb(id) = 2*tol;
                case 3 % Force backward flux
                    m.ub(id) = -2*tol;
            end
        end
    end
end