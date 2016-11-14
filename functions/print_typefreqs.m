function print_typefreqs(EEG)
% print # and % of each trial type remaining after rejects

% types of rejection
rej_type = {'rejmanual', 'rejjp', 'rejkurt', 'rejthresh'};

% get unique types
ev_types = unique({EEG.event.type}, 'stable')';
n_ev_types = length(ev_types);

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

% print out the results
fprintf('\n');
fprintf('---------------------------------'); fprintf('\n');
fprintf('event        |    # /  tot (%%)  '); fprintf('\n');
fprintf('---------------------------------'); fprintf('\n');
for ev_t = 1:n_ev_types
    type_vec = strcmpi( ep_type, ev_types{ev_t} );
    type_count = sum(type_vec);
    type_kept = sum(type_vec & ~rej_vec);
    type_percent = type_kept ./ type_count .* 100;
    type_string = sprintf('%-13s| %4d / %4d (%3.0f%%)', ev_types{ev_t}, type_kept, ...
        type_count, type_percent);
    disp(type_string);
end
fprintf('-------------------------------'); fprintf('\n');
total_rej = sum(~rej_vec);
percent_rej = total_rej ./ n_ep .* 100;
type_string = sprintf('%-13s| %4d / %4d (%3.0f%%)', 'overall', total_rej, ...
        n_ep, percent_rej);
disp(type_string);

end