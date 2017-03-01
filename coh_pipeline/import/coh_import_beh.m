function [mat_list, behdata] = coh_import_beh(opt, demogsfile)
% Imports computed behavioral measures into the MATLAB workspace for freestyle
% plotting and exploration.

% make a list of each subject's .mat file
if nargin < 2 % if no demogsfile was specified, use this
    mat_list=coh_updatemats(opt.outpath);
else
    mat_list=coh_updatemats(opt.outpath,demogsfile);
end


ns=length(mat_list);

datatype='behav_data';

%pre-allocate vars
behdata=[];

%initialize a bad subject counter(?)
bad_s=zeros(ns,1);

%load in data
s_valid=0;
for s_attempt=1:ns
load(mat_list{s_attempt},datatype);
s_valid=s_valid+1;
if exist('behav_data','var')
    behdata(s_valid)=behav_data;
end
fprintf(num2str(s_attempt))
if mod(s_attempt,20)==0
    fprintf('\n')
end
end
fprintf('\n')

if exist('behdata','var')
    beh_import
end

end