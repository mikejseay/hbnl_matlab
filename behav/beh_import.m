function behdata=beh_import(behmat_path,mat_list)
%  mobilize behavioral data into better matrices

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
    %if any(any(isnan(behav_data.respRT)))
    %    bad_s(s_attempt)=1;
    %    continue
    %end
end
s_valid=s_valid+1;
if exist('behav_data','var')
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