function coh_makesuperbatches(opt,simul_batches,version,platform)

if nargin<4
    platform='matlab';
end

if strcmpi(platform,'matlab')
    unix_bool=false;
elseif strcmpi(platform,'unix')
    unix_bool=true;
end

if nargin<3
    version=2014;
end

if version==2014
    matlabpath_prefix='/export/matlab2014a';
elseif version==2015
    matlabpath_prefix='/software/matlab2015b';
end

if isfield(opt, 'cleanset')
    if opt.cleanset
        matlabcommand=[matlabpath_prefix, '/bin/matlab -nodisplay'];
    else
        matlabcommand=[matlabpath_prefix, '/bin/matlab -nojvm -nodisplay'];
    end
else
    matlabcommand=[matlabpath_prefix, '/bin/matlab -nojvm -nodisplay'];
end
    
n_batches=length(opt.batchscripts);
n_sets=ceil(n_batches/simul_batches);

if ~unix_bool
    
    %matlab style
    for set=1:n_sets
        start_ind=(set-1)*simul_batches+1;
        if set==n_sets && set>1
            end_ind=start_ind+rem(n_batches,simul_batches)-1;
        else
            %if 
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

else
    
    %shell script style
    filename=[opt.batch_id,'_super.sh'];
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
    fprintf(fid, '#!/bin/sh\n\n');

    for command=1:n_batches
        
        if mod(command,simul_batches)==0
            suffix='';
        else
            suffix=' &';
        end
        
        [~,vlogname,~]=fileparts(opt.logpath{command});
        vlogname=[vlogname,'.vlog'];
        vlogpath=fullfile(opt.batchpath,vlogname);

        cform=[matlabcommand,' < ',opt.batchscripts{command},' > ',vlogpath,suffix];

        fprintf(fid, [cform,'\n']);
    end
    
    fclose(fid);
    
end

end