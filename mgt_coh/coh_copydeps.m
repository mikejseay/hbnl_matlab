%%

target_dir='/export/home/mike/coh_analysis';

for d=1:length(deps)
    [folder,filename,ext]=fileparts(deps{d});
    if exist(filename,'builtin')==5
        continue
    else
        copyfile(deps{d},fullfile(target_dir,[filename,ext]));
    end
end