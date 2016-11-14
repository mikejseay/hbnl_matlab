function [mat_list, id_list] = coh_updatemats(matpath, demogsfile, use_ext)
%update mat files from directory

if nargin<3
    use_ext = '.mat';
end

if nargin<2
    demogsfile=[];
end

if iscell(matpath)
    mat_list=matpath;
else

% create cell array of full paths
list=dir(matpath);
md=0;
mat_list=cell(1,1);
for file=1:length(list)
    [~,~,ext]=fileparts(list(file).name);
    if strcmpi(ext,use_ext)
        md=md+1;
        mat_list{md}=fullfile(matpath,list(file).name);
    end
end
end

%optionally reduce the list to only a certain unique file string list
if exist(demogsfile,'file')

    load(demogsfile)
    
    if exist('demogs_table','var')
        uniquefilestring=demogs_table.UniqueFileString;
    end
else % attempt to just extract the ids from the first n letters
    
    [~, id_file_list] = cellfun( @fileparts, mat_list, 'uni', 0);
    id_list = cellfun( @(x)x(1:8), id_file_list, 'uni', 0);
    return
    
end
    
n_files=length(mat_list);

[dir_string,~,~]=fileparts(mat_list{1});
string_loc=length(dir_string)+8;
%get just the unique file strings
file_string=cell(n_files,1);
for f=1:n_files
    file_string_tmp=mat_list{f}(string_loc:string_loc+10);
    file_string_tmp(2) = '1';
    file_string{f} = file_string_tmp;
end

[~,is1,~]=intersect(file_string, uniquefilestring);    
%file_string_is={file_string{is1}}';
%uniquefilestring_is={uniquefilestring{is2}}';

%mat_list_table=table(file_string_is);
%demogs_table=table(uniquefilestring_is);    

%[~,subset_inds]=join(demogs_table,mat_list_table,'LeftKeys','uniquefilestring','RightKeys','file_string');

mat_list={mat_list{is1}}';
id_list = {file_string{is1}}';
    
end