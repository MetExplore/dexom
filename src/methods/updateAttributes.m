function options = updateAttributes(prevOpt, newOpt, overwrite)
    options = prevOpt;
    fields = fieldnames(newOpt);
    for i = 1:length(fields)
        if overwrite || ~isfield(options, fields{i})
            options = setfield(options, fields{i}, getfield(newOpt, fields{i}));   
        end
    end
end