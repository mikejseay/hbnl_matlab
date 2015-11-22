%create demographic columns corresponding to the mat_list

%%

s_inds_g=false(length(s_inds),11);


%%
%load the mat containing the longer list of demographic info
%load('/export/home/mike/matlab/fmri/fmri_subs_withsession.mat')
%load('/export/home/mike/matlab/mgt_coh/ern_demog_ctl.mat')
%load('/export/home/mike/matlab/fmri/fmri_gng_subs_demogs.mat')
load('/export/home/mike/matlab/origin/fmri/fmri_demogs_simple.mat')

n_files=length(mat_list);

%discover character location of first letter of the unique file string
[dir_string,~,~]=fileparts(mat_list{1});
string_loc=length(dir_string)+8;
%get just the unique file strings
file_string=cell(n_files,1);
for f=1:n_files
    file_string{f}=mat_list{f}(string_loc:string_loc+10);
end

%%

if true

%find strings present in both
[~,is1,is2]=intersect(file_string, uniquefilestring);

%create tables and match indices with table join
%bigtable=table(SortedFileString,DoB,EEG_Age,Rec_Date,Gender);
%bigtable=table(SortedFileString,Group);
uniquefilestring=uniquefilestring(is2); group=group(is2); age_eeg=age_eeg(is2);
bigtable=table(uniquefilestring,group,age_eeg);
%bigtablename=nominal(SortedFileString);
%bigtable=dataset({SortedFileString,'ID'},{Group,'Group'});
file_string=file_string(is1);
smalltable=table(file_string);
%smalltablename=nominal(file_string);
%smalltable=dataset({file_string,'ID'});
%s_demogs=join(smalltable,bigtable);
s_demogs=join(smalltable,bigtable,'LeftKeys','file_string','RightKeys','uniquefilestring');
%s_demogs=join(smalltable,bigtable,'LeftKeys','file_string','RightKeys','SortedFileString');

end

%cut out subjects that had insufficient data
s_demogs(find(imp.bad_s==1),:)=[];

if false
%create an index matrix for age
s_median_age=median(s_demogs.EEG_Age);
s_inds_g(:,1) = (s_demogs.EEG_Age < s_median_age & s_demogs.Gender == 1) & s_inds;
s_inds_g(:,2) = (s_demogs.EEG_Age > s_median_age & s_demogs.Gender == 1) & s_inds;
s_inds_g(:,3) = (s_demogs.EEG_Age < s_median_age & s_demogs.Gender == 2) & s_inds;
s_inds_g(:,4) = (s_demogs.EEG_Age > s_median_age & s_demogs.Gender == 2) & s_inds;

s_inds_g(:,5) = s_demogs.Gender == 1 & s_inds;
s_inds_g(:,6) = s_demogs.Gender == 2 & s_inds;

s_inds_g(:,7) = s_demogs.EEG_Age < s_median_age & s_inds;
s_inds_g(:,8) = s_demogs.EEG_Age > s_median_age & s_inds;
end

%%
s_inds_g(:,9) = s_inds;

if false
s_inds_g(:,10) = s_inds & strcmpi(s_demogs.group,'Comparison');
s_inds_g(:,11) = s_inds & strcmpi(s_demogs.group,'Alcoholic');
fprintf('%d CTL subjects and %d ALC subjects remain\n',sum(s_inds_g(:,10)),sum(s_inds_g(:,11)))
end

%age_label={['Younger (< ',num2str(round(s_median_age)),')'],...
    %['Older (> ',num2str(round(s_median_age)),')']};
%age_label={['< ',num2str(round(s_median_age))],['> ',num2str(round(s_median_age))],'All'};
g_label={'Younger Male','Younger Female','Older Male','Older Female', ...
    'Male','Female','Younger','Older','All Subjects','CTL','ALC'};
g_color={'r','g','b','c','r','g','r','g','k','g','r'};
g_all=9;

nameOfStruct2Update='scl';
scl=v2struct(g_label,g_color,g_all,nameOfStruct2Update,{'fieldNames','g_label','g_color','g_all'});
clear g_label g_color nameOfStruct2Update g_all

clear n_files string_loc s_median_age DoB EEG_Age f Gender Rec_Date ...
    SortedFileString bigtable file_string Group smalltable s_inds ...
    age_eeg dir_string group uniquefilestring is1 is2