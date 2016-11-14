function exprej_struct = save_exprejects(EEG)

% types of rejection
rej_type = {'rejmanual', 'rejjp', 'rejkurt', 'rejthresh'};

% total events
n_ep = length(EEG.epoch);

% initialize logical matrix
rej_log = false(n_ep, length(rej_type));

% compile all rejects into one logical vector
rej_struct = EEG.reject;
for rej_method = 1:length(rej_type)
    if isfield(rej_struct, rej_type{rej_method}) 
        rej_contents = rej_struct.(rej_type{rej_method});
        if ~isempty(rej_contents)
            rej_log(:, rej_method) = logical(rej_contents);
        end
    end
end
rej_vec = any(rej_log, 2); % combine the rejects across types in an OR fashion

% prepare a vector of the first type code in each epoch
ep_type = cell(n_ep, 1);
for ep = 1:n_ep
    ep_type{ep} = EEG.epoch(ep).eventtype{1};
end

% find unique experiments
ep_types_first3 = cellfun( @(x)x(1:3), ep_type, 'uni', 0);
ep_exps = unique(ep_types_first3, 'stable');

% build a structure that keeps each experiment's events separate
exp_rejs = cell(length(ep_exps), 1);
exprej_struct = struct('name', ep_exps, 'rejs', exp_rejs);

for exp=1:length(ep_exps)
    tmp_evs = strcmpi(ep_types_first3, ep_exps{exp});
    exprej_struct(exp).rejs = rej_vec( tmp_evs );
end

end