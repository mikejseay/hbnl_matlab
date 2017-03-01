%% bit to plot out bad data

%2,5,17,25, and 32 are examples of high coherence values (>.999)
%3 is an example of highly correlated coherence values

for s=pp.chosen_s
figure(s)
    for chanpair=1:imp.maxpairs
        plot(cohdata(:,9,1,chanpair,s),'Color',scl.p_color(chanpair,:),'LineWidth',2); hold on
    end
hold off
legend(scl.p_label)
end

%% look at all frequencies for a given condition and pair

chosen_pair=3;

figure;
for freq=pp.chosen_freq
    plot(mean(abs(cohdata(:,freq,pp.chosen_cond,chosen_pair,pp.chosen_s)),2),'Color',scl.f_color(freq,:)); hold on
    %plot(mean(coh_results_all_test(:,freq,pp.chosen_cond,chosen_pair),2),'Color',scl.f_color(freq,:)); hold on
    %plot(mean(cohdata(:,freq,pp.chosen_cond,chosen_pair,pp.chosen_s),2),'Color',scl.f_color(freq,:)); hold on
end
hold on; plot(ones(imp.maxfreqs,1)*scl.t_zero,linspace(0,1,imp.maxfreqs),'k--'); hold off;
axis([scl.t_start scl.t_end pp.coh_absmin 1]);
set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
legend(scl.f_nlabel{pp.chosen_freq})
title(['Subject ',scl.s_label{pp.chosen_s},' / ',scl.cond_label{pp.chosen_cond},' / ',scl.p_label{chosen_pair}])

%% look at all scl.freqs for all conditions at one pair

chosen_pair=3;

figure;
for cond=1:imp.maxconds
    sp(cond)=subplot(sp_d(1),sp_d(2),cond);
    for freq=pp.chosen_freq
        plot(mean(abs(cohdata(:,freq,cond,chosen_pair,pp.chosen_s)),2),'Color',scl.f_color(freq,:)); hold on
        %plot(mean(cohstats(:,freq,cond,chosen_pair,pp.chosen_s),2),'Color',scl.f_color(freq,:)); hold on
        %plot(mean(coh_results_all_test(:,freq,cond,chosen_pair),2),'Color',scl.f_color(freq,:)); hold on
    end
    hold off; grid on; axis([scl.t_start scl.t_end pp.coh_absmin 1]);
    %hold off; grid on; axis([scl.t_start scl.t_end 0 16]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
    title(['Subject ',scl.s_label{pp.chosen_s},' / ',scl.p_label{chosen_pair},' / ',scl.cond_label{cond}])
end
linkaxes(sp)
legend(scl.f_nlabel{pp.chosen_freq})

%% look at scl.freqs for all conditions at one pair, for all S's

chosen_pair=3;

figure;
for cond=1:imp.maxconds
    sp(cond)=subplot(sp_d(1),sp_d(2),cond);
    for freq=pp.chosen_freq
        for s=1:imp.s_valid
            %plot(mean(atanh(abs(cohdata(:,freq,cond,chosen_pair,pp.chosen_s))),2),'Color',scl.s_color(pp.chosen_s,:)); hold on
            plot(mean(cohdata(:,freq,cond,chosen_pair,s),2),'Color',scl.s_color(s,:)); hold on
            %plot(mean(cohstats(:,freq,cond,chosen_pair,pp.chosen_s),2),'Color',scl.s_color(pp.chosen_s,:)); hold on
        end
    end
    %hold off; grid on; axis([scl.t_start scl.t_end .75 4]);
    hold off; grid on; axis([scl.t_start scl.t_end pp.coh_absmin 1]);
    set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
    title([scl.cond_label{cond},' / ',scl.f_nlabel{pp.chosen_freq},' Hz / ',scl.p_label{chosen_pair}])
end
linkaxes(sp)
legend(scl.s_label)

%% look at scl.freqs for all conditions at one pair, averaged across all S's

pair_subset=[4 11 14];

for freq=pp.chosen_freq
    figure;
    for cond=1:imp.maxconds+1
        sp(cond)=subplot(sp_d(1)+1,sp_d(2),cond);
        for pair=pair_subset
            if cond==imp.maxconds+1
                coh_plot_data = mean(mean(mean(cohdata(:,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),3),4),5) - ...
            mean(mean(mean(cohdata(:,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),3),4),5);
                coh_plot_data_std = std(mean(mean(cohdata(:,freq,pp.cond_diff{1},pair,s_inds_g(:,scl.g_all)),3),4) - ...
            mean(mean(cohdata(:,freq,pp.cond_diff{2},pair,s_inds_g(:,scl.g_all)),3),4),0,5); %/sqrt(sum(s_inds_g(:,scl.g_all)));
            else
                coh_plot_data=mean(mean(cohdata(:,freq,cond,pair,s_inds_g(:,scl.g_all)),4),5);
                coh_plot_data_std=std(mean(cohdata(:,freq,cond,pair,s_inds_g(:,scl.g_all)),4),0,5); %/sqrt(sum(s_inds_g(:,scl.g_all)));
            end
            plot(coh_plot_data,'Color',scl.p_color(pair,:)); hold on; %,'Color',scl.s_color(pp.chosen_s,:)); hold on
            %shadedErrorBar(1:imp.maxtimepts,coh_plot_data,coh_plot_data_std); hold on;
        end
        if cond==imp.maxconds+1
            grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit_diff(1) pp.chosen_coh_limit_diff(2)]);
        else
            grid on; axis([scl.t_start scl.t_end pp.chosen_coh_limit(1) pp.chosen_coh_limit(2)]);
        end
        plot(ones(10,1)*scl.t_zero,linspace(0,1,10),'k--'); hold off;
        set(gca,'XTick',scl.t_xtick,'XTickLabel',scl.t_xtick_ms)
        title([scl.cond_label{cond},'/',scl.f_nlabel{freq},'Hz'])
        clickableLegend(scl.p_label(pair_subset))
    end
tightfig; dragzoom;
end
distFig('s','ext');
