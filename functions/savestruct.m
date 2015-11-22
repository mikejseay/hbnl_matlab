function savestruct(outfile,instruct,version)
%function to take the contents of a struct, serialize all of its contents,
%and save them together into one mat

outstruct=[];

field_names=fieldnames(instruct);
n_fields=length(field_names);

for f=1:n_fields
    field_contents=getfield(instruct,field_names{f});
    field_contents_ser=hlp_serialize(field_contents);
    outstruct=setfield(outstruct,field_names{f},field_contents_ser);
    clear field_contents field_contents_ser
end

save(outfile,'-struct','outstruct',version)

end