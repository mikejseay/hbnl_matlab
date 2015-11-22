function coh_makesuperbatches(opt,simul_batches)

if opt.cleanset
    matlabcommand='/export/matlab2014a/bin/matlab -nodisplay';
elseif ~opt.cleanset
    matlabcommand='/export/matlab2014a/bin/matlab -nojvm -nodisplay';
end

n_batches=length(opt.batchscripts);
n_sets=ceil(n_batches/simul_batches);

for set=1:n_sets
    start_ind=(set-1)*simul_batches+1;
    if set==n_sets && set>1
        end_ind=start_ind+rem(n_batches,simul_batches)-1;
    else
        end_ind=set*simul_batches;
    end
    
    filename=[opt.batch_id,'_super',num2str(set),'.m'];
    filepath=fullfile(opt.batchpath,filename);
    
    if exist(filepath,'file')
        overwrite=input('Superbatch being overwritten. Continue? (y/n) ','s');
        if strcmpi(overwrite,'y')
            delete(filepath);
        else
            fprintf('Exiting...\n');
            return
        end
    end
    
    fid=fopen(filepath,'w');
    
    for command=start_ind:end_ind
        
        [~,vlogname,~]=fileparts(opt.logpath{command});
        vlogname=[vlogname,'.vlog'];
        vlogpath=fullfile(opt.batchpath,vlogname);
        
        cform=[matlabcommand,' < ',opt.batchscripts{command},' > ',vlogpath,' &'];
        
        fprintf(fid, ['system(''',cform,''');\n']);
    end
    
    fclose(fid);
    
end

end