function [scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile)
% index-match demographic information to eeg data matrices

% first reject subjects based on number of trials per condition
s_inds = coh_reject (imp, n_trials_all);

% create demographic columns corresponding to the mat_list
s_inds_g=false(length(s_inds),3);

load(demogsfile)

n_files=length(mat_list);

%discover character location of first letter of the unique file string
[dir_string,~,~]=fileparts(mat_list{1});
string_loc=length(dir_string)+8;
%get just the unique file strings
file_string=cell(n_files,1);
for f=1:n_files
    file_string{f}=mat_list{f}(string_loc:string_loc+10);
end

if true
%find strings present in both
[~,is1,is2]=intersect(file_string, uniquefilestring);

uniquefilestring=uniquefilestring(is2); group=group(is2); age_eeg=age_eeg(is2);
bigtable=table(uniquefilestring,group,age_eeg);
file_string=file_string(is1);
smalltable=table(file_string);
s_demogs=join(smalltable,bigtable,'LeftKeys','file_string','RightKeys','uniquefilestring');
end

%cut out subjects that had insufficient data
s_demogs(find(imp.bad_s==1),:)=[];

% start filling in indices
s_inds_g(:,1) = s_inds;

if true
s_inds_g(:,2) = s_inds & strcmpi(s_demogs.group,'Comparison');
s_inds_g(:,3) = s_inds & strcmpi(s_demogs.group,'Alcoholic');
fprintf('%d CTL subjects and %d ALC subjects remain\n',sum(s_inds_g(:,2)),sum(s_inds_g(:,3)))
end

g_label={'All Subjects','CTL','ALC'};
g_color={'b','g','r'};
g_all=1;

nameOfStruct2Update='scl';
scl=v2struct(g_label,g_color,g_all,nameOfStruct2Update,{'fieldNames','g_label','g_color','g_all'});

end