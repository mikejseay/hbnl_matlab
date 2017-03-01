function behdata_all = behav_importfolder( target_dir, expstruct, exps )
% given a subject folder, extract behavioral info of tasks in expstruct

% index its contents, pull out names only
contents = dir(target_dir);
file_names = { contents.name }';

% find .cnt files, extract their experiment names
cnt_inds = find( ~cellfun('isempty',regexp(file_names, 'cnt')) & ...
    cellfun('isempty',regexp(file_names, 'orig')) );
cnt_names = file_names( cnt_inds );
experiment_names_cnt = cellfun( @(x)x(1:3), cnt_names, 'uni', 0);

% get indices of experiment to target
if isempty(exps)
    exp_ind_set = [1:length(expstruct)]';
else
    [~, exp_ind_set] = intersect({expstruct.name}, exps);
    exp_ind_set = sort(exp_ind_set);
end

% for each experiment
edum=0;
for exp_ind = exp_ind_set'
    edum = edum + 1;
    
    % if that experiment is missing, skip it
    if ~ismember(expstruct(exp_ind).name, experiment_names_cnt)
        continue
    end
    
    cfile = find( strcmpi(experiment_names_cnt, expstruct(exp_ind).name) );
    target_cntfile = cnt_names{ cfile };
    target_cntpath = fullfile( target_dir, target_cntfile );
    
    % extract event table from CNT
    etable = cnt_getbehav( target_cntpath, expstruct(exp_ind) );
    
    % interpret behavioral data
    behdata = beh_etable( etable, expstruct(exp_ind) );
    
    % save into structure
    behdata.file = target_cntpath;
    if edum==1
        behdata_all = behdata;
    else
        behdata_all = [behdata_all, behdata];
    end

end