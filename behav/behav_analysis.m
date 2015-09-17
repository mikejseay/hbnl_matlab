%  mobilize behavioral data into better matrices

ns=length(behdata);
%n_conds=size(behdata(1).respRT,2);
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

behdataB=v2struct({'fieldNames',field_names{:}});

clearvars(field_names{:})

clear ns field_names f field_contents dims formatter d s