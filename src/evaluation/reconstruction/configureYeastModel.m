function model = configureYeastModel(yeastModel, constraintMode, minGrowthFlux)
    model = yeastModel;
    exchangeRxns = findExcRxns(model);
    % If constraintMode is provided, constrain the model to be
    % aerobic/anaerobic. Otherwise start with all the exchanges open.
    
    if constraintMode == 0
        model.lb(exchangeRxns) = -1000;
        model.ub(exchangeRxns) = 1000;
    end
    
    % Minimal medium
    if constraintMode >= 1
        model.lb(exchangeRxns) = 0;
        model.ub(exchangeRxns) = 1000;
        desiredExchanges = {...
            'r_1654'; ... % 'ammonium exchange';
            'r_1992'; ... % 'oxygen exchange'; 
            'r_2005'; ... % 'phosphate exchange';
            'r_2060'; ... % 'sulphate exchange'; 
            };
        glucoseExchange = {...
            'r_1714'; ... % D-glucose exchange'
            };
        uptakeRxnIndexes = findRxnIDs(model,desiredExchanges);
        glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
        model.lb(uptakeRxnIndexes(uptakeRxnIndexes > 0))=-1000;
        model.lb(glucoseExchangeIndex)=-10;
    end

    % Anaerobic constraints (minimal medium + sterols)
    if constraintMode == 2
        model.lb(findRxnIDs(model, 'r_1992')) = 0; % close oxygen
        sterolExchanges = {...
                'r_1757'; ... % 'ergosterol exchange';
                'r_1753'; ... % 'episterol exchange'; 
                'r_1788'; ... % 'fecosterol exchange';
                'r_1915'; ... % 'lanosterol exchange'; 
                'r_2106'; ... % 'zymosterol exchange';
                'r_2137'; ... % 'ergosta-5,7,22,24(28)-tetraen-3beta-ol exchange';
                'r_2134'; ... % '14-demethyllanosterol exchange'
                };
        uptakeSterol = findRxnIDs(model, sterolExchanges);
        model.lb(uptakeSterol) = -1000;      
    end
    
    % Force flux through the growth reaction if indicated
    model.lb(find(strcmp(model.rxnNames, 'growth'))) = minGrowthFlux;
end