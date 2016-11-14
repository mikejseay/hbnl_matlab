function [EEG, rejstruct] = save_rejects(EEG)
% save existing reject marks in an EEGLAB dataset, as a workspace variable
% and as an indepedent field of the EEG structure

% types of rejection
rej_type = {'rejmanual', 'rejjp', 'rejkurt', 'rejthresh'};

rejstruct = struct();

EEG.indreject = EEG.reject;

for rt = 1:length(rej_type)
    contents = EEG.reject.(rej_type{rt});
    rejstruct.(rej_type{rt}) = contents;
    EEG.indreject.(rej_type{rt}) = contents;
end

end