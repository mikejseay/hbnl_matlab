% check_enigma_compat_57.m - For 57 subjects, compares the enigma
% calculation results starting from two different original files: those
% that use MassComp calculations and those that have more relatively "raw"
% data from Neuroscan.

%load the 57 subjects info
load('enigma_57subs.mat')

%dummy
good_s=0;
bad_s=0;
mc_bad=0;
ns_bad=0;

clear good_ss

%for each
for s=1:length(subjectid_table)

%load MassComp data, do enigma analysis, and save output as outputvar1
mc_out=enigma_analysis_mc(mc_files{s});

%load NeuroScan data, do enigma analysis, and save output as outputvar2
ns_out=enigma_analysis_ns(ns_files{s});

%if we reached this step, it was a good subject so increment the counter
if isstruct(mc_out) && isstruct(ns_out)
    good_s=good_s+1;
    good_ss(good_s)=s;
elseif isempty(mc_out)
    mc_out=0;
    mc_bad=mc_bad+1;
elseif isempty(ns_out)
    ns_out=0;
    ns_bad=ns_bad+1;
else
    bad_s=bad_s+1;
    bad_trialnum(bad_s,1)=mc_out;
    bad_trialnum(bad_s,2)=ns_out;
end

end