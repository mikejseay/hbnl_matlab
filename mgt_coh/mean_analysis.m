%% error bar plot of a TFROI for ITC in a region
% original

t_ms = [200 400];
f_hz = [2 3];

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

[~,t_start] = min( abs(scl.t_ms - t_ms(1) ));
[~,t_end] = min( abs(scl.t_ms - t_ms(2) ));

[~,f_start] = min( abs(scl.freqs - f_hz(1) ));
[~,f_end] = min( abs(scl.freqs - f_hz(2) ));

bar_plot_y = zeros(1,2);
bar_plot_e = zeros(1,2);

for hyp=1:3
%figure;
plot_hypinds=find(opt.pair_inds==hyp);
plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));

for group=pp.plotn_g

itc_plot_base = meanx( wave_evknormdata(t_start_b:t_end_b, plot_hypchans, ...
    f_end:f_start, :, s_inds_g(:,group)), 5 );
clear itc_diff_data
for cond=1:2
    itc_diff_data{cond}=bsxfun(@minus, ...
        meanx( wave_evknormdata(t_start:t_end, plot_hypchans, f_end:f_start, ...
            cond, s_inds_g(:,group)), 5 ), ...
        itc_plot_base);
    bar_plot_y(cond,hyp) = mean( itc_diff_data{cond} );
    bar_plot_e(cond,hyp) = std( itc_diff_data{cond} ) / sqrt(sum( s_inds_g(:,group) ));
end
%errorbarbar(bar_plot_y, bar_plot_e);
%bar(bar_plot_y(hyp,:)); hold on;
%errorb(bar_plot_y(hyp,:), bar_plot_e(hyp,:));
%hold off;
end
end
%set_print_size
%tightfig;
clear_plotassistvars

%% error bar plot of a TFROI for ITC in a region
% groups

xls_dir = '/export/home/mike/matlab/xls/';
xls_filename = 'err_nki_bars.xlsx';
xls_column = 'B';
xls_row = 3;
xls_prefix = 'nki';
measure_type = 'itc';

%spec_chans = [7, 25, 58];
t_ms = [200 400];
f_hz_s = [2 4 8];
f_hz_e = [3 6 13];
f_hz_name = {'d','th','a'};

signal_type = 'region';

[~,t_start_b]=min(abs(scl.t_ms-pp.t_start_b_ms));
[~,t_end_b]=min(abs(scl.t_ms-pp.t_end_b_ms));

[~,t_start] = min( abs(scl.t_ms - t_ms(1) ));
[~,t_end] = min( abs(scl.t_ms - t_ms(2) ));

bar_plot_y = zeros(1,2);
bar_plot_e = zeros(1,2);

plotdata_cell = cell(2,1);
for freq_range=2 %:3
[~,f_start] = min( abs(scl.freqs - f_hz_s(freq_range) ) );
[~,f_end] = min( abs(scl.freqs - f_hz_e(freq_range) ) );
xls_row_count = xls_row;
for group=pp.chosen_g(pp.plotn_g)
for hyp=1:3
%plot_hypinds=find(opt.pair_inds==hyp);
%plot_hypchans=unique(opt.coherence_pairs(plot_hypinds,:));
plot_hypchans=spec_chans(hyp);

    
switch measure_type
    case 'itc'
    plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, ...
        meanx(wave_evknormdata(t_start:t_end,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)), [4 5])' ];
    
    case 'itcbl'
    plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, meanx(bsxfun(@minus, ...
        wave_evknormdata(t_start:t_end,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,plot_hypchans,f_end:f_start,:,s_inds_g(:,group)),1),4)),[4 5])' ];
    
    case 'icps'
    plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, ...
        meanx(cohdata(t_start:t_end,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)), [3 5])' ];
    
    case 'icpsbl'
    plotdata_cell{group-1, freq_range} = [ plotdata_cell{group-1, freq_range}, meanx(bsxfun(@minus, ...
        wave_evknormdata(t_start:t_end,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)), ...
        mean(mean(wave_evknormdata(t_start_b:t_end_b,f_end:f_start,:,plot_hypinds,s_inds_g(:,group)),1),3) ) ,[3 5])' ];
    
end

end
sheetname = strjoin( {xls_prefix, measure_type, f_hz_name{freq_range}, signal_type}, '_');
xls_range = [ xls_column, num2str(xls_row_count) ];
xlswrite( fullfile(xls_dir, xls_filename), plotdata_cell{group-1, freq_range}, ...
    sheetname, xls_range);
xls_row_count = xls_row_count + sum(s_inds_g(:,group));
end
end
clear_plotassistvars
