function behdata=beh_import(behmat_path,mat_list)
% import behavioral data from a folder of .mat's, each containing a structure
% called 'behav_data'

% written by michael seay, hbnl, 2015

behmat_list=coh_updatemats(behmat_path);
behmat_list=match_uniqids(behmat_list,mat_list);

%declare number of subjects to look at
ns=length(behmat_list);

%indicate mats to load
datatypes={'behav_data'};

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%load in data
s_valid=0;
for s_attempt=1:ns
    for datatype=1:length(datatypes)
        load(behmat_list{s_attempt},datatypes{datatype});
    end
    if exist('behav_data','var')
        s_valid=s_valid+1;
        behdata(s_valid)=behav_data;
    end
end
fprintf('\n')

ns=length(behdata);
field_names=fieldnames(behdata);

for f=1:length(field_names)
    
    field_contents=getfield(behdata,field_names{f});
    dims=size(squeeze(field_contents));
    formatter=':,';
    for d=1:length(dims)-1
        formatter=[formatter,':,'];
    end
    dims(end+1)=ns;
    eval([field_names{f},'=zeros(dims);'])
    
    
    for s=1:ns
        
        eval([field_names{f},'(',formatter,'s)=v2struct(behdata(s),{''fieldNames'',''',field_names{f},'''});'])
        
    end
    eval([field_names{f},'=squeeze(',field_names{f},');'])
    
end

behdata=v2struct({'fieldNames',field_names{:}});

clearvars(field_names{:})

behdata = behav_fixnames(behdata);

behdata=setfield(behdata,'cond1',{'P10','P50','N10','N50'});
behdata=setfield(behdata,'cond2',{'-80','-30','0','+20','+80'});
behdata=setfield(behdata,'resp',{'10','50'});

end