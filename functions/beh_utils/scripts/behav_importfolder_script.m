% VP3, AOD, GNG, CPT, ERN, ANT, STP
% define set of experiments to extract
exps = {'vp3','aod','gng','cpt','ern','ant','stp'};

% get expstruct that defines experiment properties
expstruct = build_expstruct;

% define folder to target
target_dir = '/raw_data/neuroscan/suny/ns205/49478002';

% extract behavioral data from cnts
behdata = behav_importfolder(target_dir, expstruct, exps);

%{
h1_names = { ...
    '/processed_data/cnt-h1-files/vp3/256/a-subjects/ns16-64/vp3_4_a1_a0001073_256_cnt.h1', ...
    '/processed_data/cnt-h1-files/aod/256/a-subjects/ns16-64/aod_5_a1_a0001073_256_cnt.h1', ...
    '/processed_data/cnt-h1-files/ern/256/a-subjects/ns16-64/ern_7_a1_a0001073_256_cnt.h1', ...
    '/processed_data/cnt-h1-files/ant/256/a-subjects/ns16-64/ant_4_a1_a0001073_256_cnt.h1'};

h1_names2 = { ...
    '/processed_data/cnt-h1-files/gng/256/a-subjects/ns16-64/gng_1_a1_a0001110_256_cnt.h1', ...
    '/processed_data/cnt-h1-files/cpt/256/a-subjects/ns16-64/cpt_2_a1_a0001110_256_cnt.h1', ...
    '/processed_data/cnt-h1-files/stp/256/a-subjects/ns16-64/stp_1_a1_a0001110_256_cnt.h1'};

target_dir = '/raw_data/neuroscan/uconn/ns113/10003051';
target_dir = '/raw_data/neuroscan/suny/ns620/49432071';
target_dir = '/raw_data/neuroscan/suny/ns164/a0001071';
target_dir = '/raw_data/neuroscan/suny/ns167/a0001073';
target_dir = '/raw_data/neuroscan/suny/ns193/a0001110';

%}