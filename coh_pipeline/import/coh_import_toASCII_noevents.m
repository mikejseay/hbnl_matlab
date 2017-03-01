function coh_import_toASCII_noevents(indir, outdir, demogsfile)

% make a list of each subject's .mat file
if nargin < 3 % if no demogsfile was specified, use this
    mat_list=coh_updatemats(indir);
else
    mat_list=coh_updatemats(indir, demogsfile);
end

ns=length(mat_list);

for s_attempt=1:ns
    load(mat_list{s_attempt}, 'dataR');
    
    [~, filename] = fileparts(mat_list{s_attempt});    
   
    name_tmp = [filename(1:20), '.txt'];
    save(fullfile(outdir, name_tmp), 'dataR', '-ascii');
end

end