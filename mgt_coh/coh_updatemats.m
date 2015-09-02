function mat_list=coh_updatemats(matpath,demogsfile)

%update mat files from directory


%filepath='/active_projects/mike/fmri_phase4_err_cleanCSD/';
if iscell(matpath)
    mat_list=matpath;
else
list=dir(matpath);
md=0;
mat_list=cell(2,1);
for file=1:length(list)
    [~,~,ext]=fileparts(list(file).name);
    if strcmpi(ext,'.mat')
        md=md+1;
        mat_list{md}=[matpath,list(file).name];
    end
end
end

%optionally reduce the list to only a certain unique file string list
if true

    load(demogsfile)
    n_files=length(mat_list);

    [dir_string,~,~]=fileparts(mat_list{1});
    string_loc=length(dir_string)+8;
    %get just the unique file strings
    file_string=cell(n_files,1);
    for f=1:n_files
        file_string{f}=mat_list{f}(string_loc:string_loc+10);
    end
    
    [~,is1,is2]=intersect(file_string, uniquefilestring);    
    %file_string_is={file_string{is1}}';
    %uniquefilestring_is={uniquefilestring{is2}}';
    
    %mat_list_table=table(file_string_is);
    %demogs_table=table(uniquefilestring_is);    
    
    %[~,subset_inds]=join(demogs_table,mat_list_table,'LeftKeys','uniquefilestring','RightKeys','file_string');
    
    mat_list={mat_list{is1}}';
    
end