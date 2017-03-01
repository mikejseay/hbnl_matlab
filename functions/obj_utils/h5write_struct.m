function h5write_struct(output_file,name,instruct)
%save all of the contents of a structure to an h5 file

field_names=fieldnames(instruct);
n_fields=length(field_names);

h5create(output_file,name,Inf,'Datatype','double');

for f=1:n_fields
    field_contents=getfield(instruct,field_names{f});
    h5write(output_file,name,field_contents);
    clear field_contents
end

end