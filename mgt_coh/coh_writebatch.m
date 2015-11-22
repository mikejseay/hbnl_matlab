function scriptpath = coh_writebatch(batchname,batchpath,optpath,batch_ind,start_ind,end_ind)

scriptname = [batchname,'_b',num2str(batch_ind),'.m'];
scriptpath = fullfile(batchpath,scriptname);

if exist(scriptpath,'file')
    overwrite=input('Batch being overwritten. Continue? (y/n) ','s');
    if strcmpi(overwrite,'y')
        delete(scriptpath);
    else
        fprintf('Exiting...\n');
        return
    end
end

fid=fopen(scriptpath, 'w');

fprintf(fid, ['load(''',optpath,''')\n']);
fprintf(fid,'\n');
fprintf(fid,'file_list=opt.infile_list;\n');
fprintf(fid,'output_dir = opt.outpath;\n');
fprintf(fid,'data_file_type=opt.data_type;\n');
fprintf(fid,'suffix=opt.outsuffix;\n');
fprintf(fid,['log_file=opt.logpath{',num2str(batch_ind),'};\n']);
fprintf(fid,'\n');
fprintf(fid,'fid = fopen(log_file, ''w'');\n');
fprintf(fid,'filecell = list2cell(file_list);\n');
fprintf(fid,'n_files = length(filecell);\n');
fprintf(fid,['for file = ',num2str(start_ind),':',num2str(end_ind),'\n']);
fprintf(fid,'    data_file = filecell{file};\n');
fprintf(fid,'    if ~exist(data_file, ''file'')\n');
fprintf(fid,'        fprintf(fid, ''not found: %%s\\n'', data_file);\n');
fprintf(fid,'        continue\n');
fprintf(fid,'    end\n');
fprintf(fid,'    base_file = basename(data_file);\n');
fprintf(fid,'    output_file = sprintf(''%%s/%%s.%%s.mat'', output_dir, base_file, suffix);\n');
fprintf(fid,'    if exist(output_file, ''file'')\n');
fprintf(fid,'        fprintf(fid, ''skipping %%s\\n'', output_file);\n');
fprintf(fid,'        continue;\n');
fprintf(fid,'    else\n');
fprintf(fid,'        fprintf(fid, ''starting %%s '',base_file );\n');
fprintf(fid,'    end\n');
fprintf(fid,'    try\n');
fprintf(fid,'        output = coh_calc(data_file, data_file_type, opt);\n');
fprintf(fid,'    catch err\n');
fprintf(fid,'        fprintf(fid, ''%%s: '', err.message);\n');
fprintf(fid,'        output = [];\n');
fprintf(fid,'    end\n');
fprintf(fid,'    if ~isempty(output)\n');
fprintf(fid,'        n_trials = sum(output.trials);\n');
fprintf(fid,'        fprintf(fid, ''%%d '', n_trials);\n');
fprintf(fid,'        fprintf(fid, ''\\n'');\n');
fprintf(fid,'        save(output_file, ''-struct'', ''output'', ''-v6'');\n');
fprintf(fid,'    else\n');
fprintf(fid,'        fprintf(fid, ''error; not saved\\n'');\n');
fprintf(fid,'    end\n');
fprintf(fid,'end\n');
fprintf(fid,'if fid > 1\n');
fprintf(fid,'    fclose(fid);\n');
fprintf(fid,'end\n');

if fid > 1
    fclose(fid);
end

end