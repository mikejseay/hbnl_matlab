% check to see how many (if any) of the h1_structs in an opt have a shorter
% continuous data than the trial mat says it should

epoch_start = -500;

%load('/export/home/mike/matlab/batch/nki_err_nonCSD/nki_err_cnth1snonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/nki_ern_nonCSD/nki_ern_cnth1snonCSD_opt.mat')
%load('/export/home/mike/matlab/batch/nki_ans_nonCSD/nki_ans_cnth1s_nonCSD_opt.mat')
load('/export/home/mike/matlab/batch/nki_stp_nonCSD/nki_stp_cnth1s_nonCSD_opt.mat')

file_list=opt.infile_list;
filecell = list2cell(file_list);
n_files = length(filecell);

bad_bool = false(n_files, 1);

for file=1:n_files
    h1_struct = read_hdf1_dataV7(filecell{file});
    rate = h1_struct.experiment_struct.rate;
    trial_index = h1_struct.trial_struct.time_offset * rate;
    if trial_index(1) - 384 < 1 %round(-epoch_start./1000*rate)
        bad_bool(file) = true;
    end
end

disp('Done!!!');