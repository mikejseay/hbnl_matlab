

%cnth1_file='/processed_data/cnt-h1-files/err/256/a-subjects/ns32-64/err_8_a1_a0001363_256_cnt.h1';
%avgh1_file='/processed_data/avg-h1-files/err/l16-h003-t75-b200/a-subjects/ns32-64/err_8_a1_a0001363_avg.h1';
%check_sub=1;
%avgh1_file='/processed_data/avg-h1-files/err/l16-h003-t75-b200/a-subjects/ns32-64/err_8_a1_a0001412_avg.h1';
%check_sub=4;
%avgh1_file='/processed_data/avg-h1-files/err/l16-h003-t75-b200/a-subjects/ns32-64/err_8_a1_a0001452_avg.h1';
%check_sub=8;
avgh1_file='/processed_data/avg-h1-files/err/l16-h003-t75-b200/a-subjects/ns32-64/err_8_a1_a0001379_avg.h1';
check_sub=8;
avgh1_struct = read_hdf1_dataV7(avgh1_file);

%%

check_cond=1;
check_chan=7;

prestim_pts=time_zero - (avgh1_struct.experiment_struct.pre_stim_time_ms / 1000 * ...
    avgh1_struct.experiment_struct.rate) + 1;

h1_erp_data=squeeze(mean_data_all(prestim_pts:end,check_chan,check_cond,check_sub));
avg_data=squeeze(avgh1_struct.data_struct.hdf1_avg_data(check_cond,check_chan,1:length(h1_erp_data)));

avg_scalefactor=mean(h1_erp_data./avg_data);
avg_scalefactor
div_factor=avgh1_struct.experiment_struct.uv_per_unit(1);
mult_factor=1/div_factor;

figure; plot(h1_erp_data,'r'); hold on; plot(avg_data*mult_factor,'g');
legend({'processed from cnt.h1',['from avg.h1 * ',num2str(round(mult_factor))]})