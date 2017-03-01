%% implement a RMANOVA to the ERP peaks/latencies, MATLAB built-ins

%input data must be s_inds x conds
%ranova_input=squeeze(mean(mean_data_all(221,7,:,s_inds),1))';
%ranova_input=peakmat(:,:,1);
%ranova_input=squeeze(mean(wavelet_tot_all(206:257,56,19,:,s_inds),1))';
ranova_input=squeeze(mean(mean(cohdata(180:231,14:16,:,4,:),1),2))';
ranova_tbl=s_demogs(:,2:end);
for cond=1:imp.maxconds
    tbl_col=table(ranova_input(:,cond),'VariableNames',{['c',num2str(cond)]});
    ranova_tbl=[ranova_tbl,tbl_col];
end
ranova_tbl(~s_inds_g(:,scl.g_all),:)=[];
preds=table([1:imp.maxconds]','VariableNames',{'Conditions'});
rm=fitrm(ranova_tbl,['c1-c',num2str(imp.maxconds),'~group+age_eeg'],'WithinDesign',preds);
ranova_results=ranova(rm);
anova_results=anova(rm);

%% RMANOVA with anova_rm

%input data must be s_inds x conds
%ranova_input=squeeze(mean(mean_data_all(221,7,:,s_inds),1))';
%ranova_input=peakmat(:,:,1);
%ranova_input=squeeze(mean(itc_results_all(206:257,56,19,:,s_inds),1))';
ranova_input=squeeze(mean(mean(cohdata(180:231,14:16,:,4,:),1),2))';
gdum=0;
for group=pp.chosen_g
    gdum=gdum+1;
    ranova_data{gdum}=ranova_input(s_inds_g(:,group),:);
end
[p,anova_table]=anova_rm(ranova_data);