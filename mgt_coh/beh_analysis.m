%% fix certain things

%behav_fixnames
behdata.crit=behdata.crit';
behdata.avgbet=behdata.avgbet';

%% plot the average bet following a given outcome

for cond=1:4
figure(1);
subplot(2,2,cond);
hist(behdata.avgbet_po(cond,:));
title(['Average Bet Following ',behdata.cond1{cond}]);
axis([10 50 0 15])
end

%% scatter the decision criterions in different ways

scatterin=behdata.crit_po;
figure(2)

subplot(2,2,1);
scatter(scatterin(3,:),scatterin(1,:))
xlabel('Small Loss'); ylabel('Small Gain');

subplot(2,2,2);
scatter(scatterin(4,:),scatterin(2,:))
xlabel('Big Loss'); ylabel('Big Gain');

subplot(2,2,3);
scatter(scatterin(1,:),scatterin(2,:))
xlabel('Small Gain'); ylabel('Big Gain');

subplot(2,2,4);
scatter(scatterin(3,:),scatterin(4,:))
xlabel('Small Loss'); ylabel('Big Loss');

%test directly for length of trains of repeats of one response
plottitle('Criterion (+ is betting 50)')

%%
figure(3)
subplot(2,1,1)
scatter(behdata.crit,mean(behdata.avgbet_po([1 3],:),1))
xlabel('Criterion'); ylabel('Avg. Bet Following Small Outcome');

subplot(2,1,2)
scatter(behdata.crit,mean(behdata.avgbet_po([2 4],:),1))
xlabel('Criterion'); ylabel('Avg. Bet Following Large Outcome');

%%
figure
subplot(2,1,1)
scatter(behdata.crit,mean(behdata.avgbet_po([1 2],:),1))
xlabel('Criterion'); ylabel('Avg. Bet Following Gain');

subplot(2,1,2)
scatter(behdata.crit,mean(behdata.avgbet_po([3 4],:),1))
xlabel('Criterion'); ylabel('Avg. Bet Following Loss');

%%

%[h,p,ci,stats]=ttest(crit(1:30),crit(31:60));

%%

figure(4)
clust_colors={'r','g','b','k'};
for cluster=1:size(output,1)
subplot(2,2,1);
scatter(behdata.avgprev_po(logical(output(cluster,:)),3),behdata.avgprev_po(logical(output(cluster,:)),1),clust_colors{cluster})
hold on
xlabel('Avg. Bet After Small Loss'); ylabel('Avg. Bet After Small Gain');

subplot(2,2,2);
scatter(behdata.avgprev_po(logical(output(cluster,:)),4),behdata.avgprev_po(logical(output(cluster,:)),2),clust_colors{cluster})
hold on
xlabel('Avg. Bet After Big Loss'); ylabel('Avg. Bet After Big Gain');

subplot(2,2,3);
scatter(behdata.avgprev_po(logical(output(cluster,:)),1),behdata.avgprev_po(logical(output(cluster,:)),2),clust_colors{cluster})
hold on
xlabel('Avg. Bet After Small Gain'); ylabel('Avg. Bet After Big Gain');

subplot(2,2,4);
scatter(behdata.avgprev_po(logical(output(cluster,:)),3),behdata.avgprev_po(logical(output(cluster,:)),4),clust_colors{cluster})
xlabel('Avg. Bet After Small Loss'); ylabel('Avg. Bet After Big Loss');
hold on
end

%% turn the "average bet following previous outcome" mat into a "proportion mat"

behdata.prop_po=(behdata.avgbet_po-10)/40;

%% RMANOVA with anova_rm

s_inds_g(:,1)=strcmp(group,'Alcoholic');
s_inds_g(:,2)=strcmp(group,'Comparison');

%ranova_input=valencemeans;
%ranova_input=avgprev2;
gdum=0;
for g=1:2
    gdum=gdum+1;
    ranova_data{gdum}=ranova_input(s_inds_g(:,g),:);
end
[p,anova_table]=anova_rm(ranova_data);

%% 2way RMANOVA with rm_anova2

%ranova_input=behdata.avgprev_po;   %columns are 10,50,-10,-50 {[1 1],[2 1],[1 2],[2 2]}
ranova_input=crit_prevoutcome;
Y=reshape(ranova_input,240,1);
S=[1:60]';
S=[S;S;S;S];
G=[ones(30,1);2*ones(30,1)];
G=[G;G;G;G];
F1=[ones(60,1);2*ones(60,1)];
F1=[F1;F1];
F2=[ones(120,1);2*ones(120,1)];
FACTNAMES={'amount','valence'};
stats=rm_anova2(Y,S,F1,F2,FACTNAMES);

%% build  a dataset for export to csv (for R, e.g)

ds_out=[Y,S,G,F1,F2];
ds_out=array2table(ds_out,'VariableNames',{'Criterion','Subject','Group','PrevAmount','PrevValence'});

%% RMANOVA with matlab built-ins

ranova_input=avgprev2;
ranova_tbl=table(group,age_eeg);
for cond=1:size(ranova_input,2)
    tbl_col=table(ranova_input(:,cond),'VariableNames',{['c',num2str(cond)]});
    ranova_tbl=[ranova_tbl,tbl_col];
end
%preds=table([1 2 1 2]',[1 1 2 2]',...
%    'VariableNames',{'amount','valence'});
preds=[-80 -30 0 30 80]';
%rm=fitrm(ranova_tbl,['c1-c',num2str(size(ranova_input,2)),'~group'],'WithinDesign',preds,'WithinModel','amount*valence');
rm=fitrm(ranova_tbl,['c1-c',num2str(size(ranova_input,2)),'~group'],'WithinDesign',preds,'WithinModel','orthogonalcontrasts');
ranova_results=ranova(rm);
anova_results=anova(rm);