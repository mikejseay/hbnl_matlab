function outstruct=transposefields(instruct)

outstruct=[];

field_names=fieldnames(instruct);
n_fields=length(field_names);

for f=1:n_fields
    field_contents=getfield(instruct,field_names{f});
    field_contents=field_contents';
    outstruct=setfield(outstruct,field_names{f},field_contents);
end

end