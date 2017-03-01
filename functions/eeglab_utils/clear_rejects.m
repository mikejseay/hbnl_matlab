function EEG = clear_rejects(EEG)
% clear existing reject marks in an EEGLAB dataset

% types of rejection
rej_type = {'rejmanual', 'rejjpE', 'rejjp', 'rejkurtE', 'rejkurt', ...
    'rejthreshE', 'rejthresh'};

for rt = 1:length(rej_type)
    if ~isempty( EEG.reject.(rej_type{rt}) )
        EEG.reject.(rej_type{rt}) = [];
    end
end

end