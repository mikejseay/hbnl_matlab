% load the opt
%load('/export/home/mike/matlab/batch/ph4_tb/ph4-baseline-dirs-gawked-dirsonlyph4_tb_opt.mat')
load('/export/home/mike/matlab/batch/ph4bl_tb/ph4-bl-rawcntfoldersph4bl_tb_opt.mat')

% import beh data as a table
beh_table = beh_import2( opt );

%%

out_name = [opt.outsuffix, '.csv'];
out_path = fullfile( opt.batchpath, out_name );
writetable( beh_table, out_path);