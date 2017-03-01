%% copy all dependencies into a target folder

target_dir='/export/home/mike/distributions/coh_pipeline';
old_path='/export/home/mike/matlab';

for d=1:length(deps)
    [folder,filename,ext]=fileparts(deps{d});
    if exist(filename,'builtin')==5 %check if built-in, this doesn't work?
        continue
    elseif ~isempty(regexp(filename,'join')) %check a specific name
        continue
    end
    %pull out the subfolders
    %[~,loc]=regexp(folder,old_path);
    pathloc=length(old_path)+1;
    newpaths=folder(pathloc:end);
    %check for the existence of the subfolders
    if ~exist(fullfile(target_dir,newpaths),'dir')
        %recursively create directory
        mkdir(fullfile(target_dir,newpaths));
    end
    copyfile(deps{d},fullfile(target_dir,newpaths,[filename,ext]));
end