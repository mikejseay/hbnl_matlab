function coh_import_toASCII(opt, outdir, conds, demogsfile)

% make a list of each subject's .mat file
if nargin < 4 % if no demogsfile was specified, use this
    mat_list=coh_updatemats(opt.outpath);
else
    mat_list=coh_updatemats(opt.outpath,demogsfile);
end

if nargin < 3
    conds=1:length(opt.case_vec);
end

ns=length(mat_list);

for s_attempt=1:ns
    load(mat_list{s_attempt}, 'erptrial', 'erptrial_log');
    
    [~, filename] = fileparts(mat_list{s_attempt});
    n_trials = sum(erptrial_log);
    
    % here you would transform it like this
    for cond=conds
        erptrial_tmp = erptrial(:, erptrial_log(:, cond), :);
        erptrial2d = reshape(erptrial_tmp, ...
            size(erptrial_tmp,1)*size(erptrial_tmp,2), ...
            size(erptrial_tmp,3));
        name_tmp = [filename(1:22), num2str(size(erptrial_tmp,1)), '_', num2str(n_trials(cond)), ...
            '_', opt.case_label{cond}, '.txt'];
        save(fullfile(outdir, name_tmp), 'erptrial2d', '-ascii');
    end
end

end