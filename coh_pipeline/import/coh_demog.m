function [scl, s_inds_g, s_demogs] = coh_demog(mat_list, imp, scl, n_trials_all, demogsfile, trials_necessary, age_range, custom_rej)
% index-match demographic information to eeg data matrices

if nargin<8
    custom_rej=[];
end
if nargin<7
    age_range=[];
end
if nargin<6
    trials_necessary=15;
end
if nargin<5 || isempty(demogsfile)
    demogs_logic=false;
    s_demogs=[];
else
    demogs_logic=true;
end

% first reject subjects based on number of trials per condition
s_inds = coh_reject (imp, n_trials_all, trials_necessary, custom_rej);

% create demographic columns corresponding to the mat_list
s_inds_g=false(length(s_inds),3);

if demogs_logic
    load(demogsfile)
    if exist('demogs_table','var')
        uniquefilestring=demogs_table.UniqueFileString;
        % group=demogs_table.POP;
        group=demogs_table.sex;
        age_eeg=demogs_table.EEG_Age;
    end
end

n_files=length(mat_list);

%discover character location of first letter of the unique file string
[dir_string,~,~]=fileparts(mat_list{1});
string_loc=length(dir_string)+8;
%get just the unique file strings
file_string=cell(n_files,1);
for f=1:n_files
    file_string_tmp=mat_list{f}(string_loc:string_loc+10);
    file_string_tmp(2) = '1';
    file_string{f}=file_string_tmp;
end

if demogs_logic
    %find strings present in both
    [~,is1,is2]=intersect(file_string, uniquefilestring);

    uniquefilestring=uniquefilestring(is2); group=group(is2); age_eeg=age_eeg(is2);
    bigtable=table(uniquefilestring,group,age_eeg);
    file_string=file_string(is1);
    smalltable=table(file_string);
    s_demogs=join(smalltable,bigtable,'LeftKeys','file_string','RightKeys','uniquefilestring');

    %cut out subjects that had insufficient data
    s_demogs(find(imp.bad_s==1),:)=[];
end

% start filling in indices
s_inds_g(:,1) = s_inds;

if demogs_logic
    
    %groups
    if iscellstr(s_demogs.group(1))
        %s_inds_g(:,2) = s_inds & strcmpi(s_demogs.group,'Comparison');
        %s_inds_g(:,3) = s_inds & strcmpi(s_demogs.group,'Alcoholic');
        %s_inds_g(:,2) = s_inds & strcmpi(s_demogs.group,'C');
        %s_inds_g(:,3) = s_inds & strcmpi(s_demogs.group,'A');
        s_inds_g(:,2) = s_inds & strcmpi(s_demogs.group,'m');
        s_inds_g(:,3) = s_inds & strcmpi(s_demogs.group,'f');
    elseif isnumeric(s_demogs.group(1))
        s_inds_g(:,2) = s_inds & s_demogs.group==1;
        s_inds_g(:,3) = s_inds & s_demogs.group==2;
    end
    
    %age mask
    if ~isempty(age_range)
        s_inds_g = bsxfun(@and, s_inds_g, s_demogs.age_eeg > age_range(1));
        s_inds_g = bsxfun(@and, s_inds_g, s_demogs.age_eeg < age_range(2));
    end
    
    fprintf('%d G1 subjects and %d G2 subjects remain\n',sum(s_inds_g(:,2)),sum(s_inds_g(:,3)))
    
end

fprintf('%d total subjects remain\n', sum(s_inds_g(:,1)) )

%g_label={'All Subjects','CTL','ALC'};
g_label={'All Subjects','Male','Female'};
g_color={'k','g','r'};
g_all=1;

nameOfStruct2Update='scl';
scl=v2struct(g_label,g_color,g_all,nameOfStruct2Update,{'fieldNames','g_label','g_color','g_all'});

end