function opt = coh_multiopt( opts )
% load multiple opts as a single opt

load( opts{1} )
tmp_opt = opt;

n_opts = length(opts);
for o = 2:n_opts
    load( opts{o} )
    old_fields = fieldnames(tmp_opt);
    new_fields = fieldnames(opt);
    if length(new_fields) > length(old_fields)
        add_fields = setdiff(new_fields, old_fields);
        for f=1:length(add_fields)
            tmp_opt(1).(add_fields{f}) = [];
        end
    else
        add_fields = setdiff(old_fields, new_fields);
        for f=1:length(add_fields)
            opt.(add_fields{f}) = [];
        end
    end
    tmp_opt = [tmp_opt, opt];
end

opt = tmp_opt;

end