function [p]=rmanova_inline(ranova_input)

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

p(1)=ranova_results(1,5);
p(2)=anova_results(2,7);
p(3)=ranova_results(2,5);

end