%% run pre-processing to get trial data, record important checkpoints

% goal is to test pre-processing methods on a small number of subjects

%load h1 file list
load('/export/home/mike/matlab/origin/mgt_coh/fmri_ern_h1_files.mat')
%params with all new steps
%load('/export/home/mike/matlab/batch/ern_fmri_pcheck.mat')
load('/export/home/mike/matlab/batch/ern_fmri_pchecknew.mat')

for i=1:length(h1_list)
[Y,name]=erp_analysis_pcheck(h1_list{i},'hdf_binary',[],param_struct,...
    [false false false false],false);
out_file=[param_struct.raw_path,'/',name,'_cleanraw.mat'];
save(out_file,'-struct','Y','-v6');
end