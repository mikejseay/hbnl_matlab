function expstruct = build_expstruct
% build struct describing experiment info

min_resptime = 100; % ms

exp_list = {'eeo', 'eec', 'vp3', 'cpt', 'ern', 'ant', 'aod', 'ans', 'stp', 'gng', 'pfb', 'err'}';
iti_list = {  0     0      1633   1650   1500   1600   1500   2000   2300   1200   2100,  1500}';

exp2ind = containers.Map(exp_list, 1:length(exp_list));

ttl_map = cell(length(exp_list), 2);
resp = cell(length(exp_list), 1);
respwin = cell(length(exp_list), 1);
init_trial = cell(length(exp_list), 1);

% eeo
ttl_map{1,1} = []; ttl_map{1,2} = {};
resp{1} = []; respwin{1} = [];
init_trial{1} = [];

% eec
ttl_map{2,1} = []; ttl_map{2,2} = {};
resp{2} = []; respwin{2} = [];
init_trial{2} = [];

% vp3
ttl_map{3,1}    = [10 20 30];
ttl_map{3,2}    = {'vp3_target', 'vp3_nt', 'vp3_novel'};
resp{3}         = [1 0 0];
respwin{3}      = [800 0 0];
init_trial{3}   = [10 20 30];

% cpt
ttl_map{4,1}    = [10 20 30 40 50 60];
ttl_map{4,2}    = {'cpt_go', 'cpt_cue', 'cpt_cuednogo', 'cpt_db4nogo', ...
                    'cpt_nogo', 'cpt_dafterd'};
resp{4}         = [1 0 0 0 0 0];
respwin{4}      = [1000 0 0 0 0 0];
init_trial{4}   = [10 20 30 40 50 60];

% ern
ttl_map{5,1}    = [80 60 20 40 90 15 75 85 95];
ttl_map{5,2}    = {'ern_n50', 'ern_n10', 'ern_p10', 'ern_p50', 'ern_choice', ...
                    'ern_noresp', 'ern_lostmsg', 'ern_evenmsg', 'ern_winmsg'};
resp{5}         = [-1 -1 -1 -1 1 -1 -1 -1 -1];
respwin{5}      = [0 0 0 0 1000 0 0 0 0];
init_trial{5}   = 90;

% ant
ttl_map{6,1}    = [10 20 30 40];
ttl_map{6,2}    = {'ant_jumble', 'ant_prime', 'ant_ant', 'ant_word'};
resp{6}         = [8 1 1 1];
respwin{6}      = [800 800 800 800];
init_trial{6}   = [10 20 30 40];

% aod
ttl_map{7,1}    = [10 20];
ttl_map{7,2}    = {'aod_target', 'aod_nt'};
resp{7}         = [1 0];
respwin{7}      = [800 0];
init_trial{7}   = [10 20];

% ans
ttl_map{8,1}    = [10 20 30 40 50 60 70 80 90];
ttl_map{8,2}    = {'ans_t1', 'ans_t2', 'ans_t3', 'ans_t4', 'ans_t5', 'ans_t6', ...
                    'ans_t7', 'ans_t8', 'ans_t9p'};
resp{8}         = [0 0 0 0 0 0 0 0 0];
respwin{8}      = [0 0 0 0 0 0 0 0 0];
init_trial{8}   = 10;

% stp
ttl_map{9,1}    = [10 20];
ttl_map{9,2}    = {'stp_cong', 'stp_incong'};
resp{9}         = [1 8];
respwin{9}      = [1000 1000];
init_trial{9}   = [10 20];

% gng
ttl_map{10,1}   = [10 20 40 50 60 70];
ttl_map{10,2}   = {'gng_go', 'gng_nogo', 'gng_dg', 'gng_xg', 'gng_xng', 'gng_dng'};
resp{10}        = [1 0 0 0 0 0];
respwin{10}     = [500 0 0 0 0 0];
init_trial{10}  = [10 20];

% pfb
ttl_map{11,1}   = [10 20 30 40 50 15];
ttl_map{11,2}   = {'pfb_lbc', 'pfb_mbc', 'pfb_rbc', 'pfb_correct', 'pfb_incorrect', ...
                        'pfb_noresp'};
resp{11}        = [1 1 1 -1 -1 -1];
respwin{11}     = [1000 1000 1000 0 0 0];
init_trial{11}  = [10 20 30];

% err
ttl_map{12,1}    = [80 60 20 40 90 15 75 85 95];
ttl_map{12,2}    = {'ern_n50', 'ern_n10', 'ern_p10', 'ern_p50', 'ern_choice', ...
                    'ern_noresp', 'ern_lostmsg', 'ern_evenmsg', 'ern_winmsg'};
resp{12}         = [-1 -1 -1 -1 1 -1 -1 -1 -1];
respwin{12}      = [0 0 0 0 1000 0 0 0 0];
init_trial{12}   = 90;

% ttl2name = cell(length(exp_list),1);
% for exp=3:length(exp_list)
%     ttl2name{exp} = containers.Map(ttl_map{exp,1}, ttl_map{exp,2});
% end

min_resptime = ones(length(exp_list), 1) .* min_resptime;

expstruct = struct('name', exp_list, 'iti', iti_list, ...
    'ttl_codes',ttl_map(:,1), 'ttl_descriptors', ttl_map(:,2), ...
    'correct_resps', resp, ...
    'resp_windows', respwin, 'min_resptime', min_resptime, ...
    'init_trial', init_trial);

end