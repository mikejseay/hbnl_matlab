function batchscripts=coh_makebatches(batch_id,infile_list,batchpath,optpath)

infile_cell=list2cell(infile_list);
n_files=length(infile_cell);

if n_files >= 400
    files_perscript = 100;
else
    files_perscript = ceil(n_files / 4);
end

n_batches=ceil(n_files/files_perscript);

batchscripts=cell(n_batches,1);

for batch=1:n_batches
    
    startfile_ind = (batch-1)*files_perscript+1;
    
    if batch==n_batches && batch > 1
        endfile_ind = startfile_ind + rem(n_files,files_perscript) - 1;
    else
        endfile_ind = batch * files_perscript;
    end
    
    batchscripts{batch}=coh_writebatch(batch_id,batchpath,optpath,batch,startfile_ind,endfile_ind);
    
end

end