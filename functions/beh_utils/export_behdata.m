% reformat a behdata structure into a table and export it as a csv

behdata.prop=(behdata.avgbet-10)/40;
behdata.prop_po=(behdata.avgbet_po-10)/40;
behdata.prop_po2=(behdata.avgbet_po2-10)/40;

file_string = s_demogs.file_string;
ID = cellfun(@(x)x(4:11), file_string, 'uni', 0);
session = cellfun(@(x)x(1), file_string, 'uni', 0);

%
choice_criterion = behdata.crit;
avgbet = behdata.avgbet;

% choice of 10, 50
bets_10 = behdata.resp_mat(1, :)';
bets_50 = behdata.resp_mat(2, :)';

% mean, median, variance
% 'Bet 10 Post-Loss','Bet 10 Post-Win','Bet 50 Post-Loss','Bet 50 Post-Win'
medianrt_bet10postloss = squeeze(behdata.respRT(2, 1, :));
medianrt_bet10postgain = squeeze(behdata.respRT(2, 2, :));
medianrt_bet50postloss = squeeze(behdata.respRT(2, 3, :));
medianrt_bet50postgain = squeeze(behdata.respRT(2, 4, :));

varrt_bet10postloss = squeeze(behdata.respRT(3, 1, :));
varrt_bet10postgain = squeeze(behdata.respRT(3, 2, :));
varrt_bet50postloss = squeeze(behdata.respRT(3, 3, :));
varrt_bet50postgain = squeeze(behdata.respRT(3, 4, :));

avgbet_postP10 = behdata.avgbet_po(1, :)';
avgbet_postP50 = behdata.avgbet_po(2, :)';
avgbet_postN10 = behdata.avgbet_po(3, :)';
avgbet_postN50 = behdata.avgbet_po(4, :)';

stdbet_postP10 = behdata.stdbet_po(1, :)';
stdbet_postP50 = behdata.stdbet_po(2, :)';
stdbet_postN10 = behdata.stdbet_po(3, :)';
stdbet_postN50 = behdata.stdbet_po(4, :)';

crit_postP10 = behdata.crit_po(1, :)';
crit_postP50 = behdata.crit_po(2, :)';
crit_postN10 = behdata.crit_po(3, :)';
crit_postN50 = behdata.crit_po(4, :)';

% -80, -30, 0, +30, +80
bets_postN80 = behdata.po_mat(1, :)';
bets_postN30 = behdata.po_mat(2, :)';
bets_postE   = behdata.po_mat(3, :)';
bets_postG30 = behdata.po_mat(4, :)';
bets_postG80 = behdata.po_mat(5, :)';

avgbet_postN80 = behdata.avgbet_po2(1, :)';
avgbet_postN30 = behdata.avgbet_po2(2, :)';
avgbet_postE   = behdata.avgbet_po2(3, :)';
avgbet_postG30 = behdata.avgbet_po2(4, :)';
avgbet_postG80 = behdata.avgbet_po2(5, :)';

stdbet_postN80 = behdata.stdbet_po2(1, :)';
stdbet_postN30 = behdata.stdbet_po2(2, :)';
stdbet_postE   = behdata.stdbet_po2(3, :)';
stdbet_postG30 = behdata.stdbet_po2(4, :)';
stdbet_postG80 = behdata.stdbet_po2(5, :)';

prop = behdata.prop;

prop_postP10 = behdata.prop_po(1, :)';
prop_postP50 = behdata.prop_po(2, :)';
prop_postN10 = behdata.prop_po(3, :)';
prop_postN50 = behdata.prop_po(4, :)';

prop_postN80 = behdata.prop_po2(1, :)';
prop_postN30 = behdata.prop_po2(2, :)';
prop_postE   = behdata.prop_po2(3, :)';
prop_postG30 = behdata.prop_po2(4, :)';
prop_postG80 = behdata.prop_po2(5, :)';


%%

beh_table = table(ID, session, choice_criterion, avgbet, bets_10, bets_50, ...
    medianrt_bet10postloss, medianrt_bet10postgain, medianrt_bet50postloss, medianrt_bet50postgain, ...
    varrt_bet10postloss, varrt_bet10postgain, varrt_bet50postloss, varrt_bet50postgain, ...
    avgbet_postP10, avgbet_postP50, avgbet_postN10, avgbet_postN50, ...
    stdbet_postP10, stdbet_postP50, stdbet_postN10, stdbet_postN50, ...
    crit_postP10, crit_postP50, crit_postN10, crit_postN50, ...
    bets_postN80, bets_postN30, bets_postE, bets_postG30, bets_postG80, ...
    avgbet_postN80, avgbet_postN30, avgbet_postE, avgbet_postG30, avgbet_postG80, ...
    stdbet_postN80, stdbet_postN30, stdbet_postE, stdbet_postG30, stdbet_postG80, ...
    prop, prop_postP10, prop_postP50, prop_postN10, prop_postN50, ...
    prop_postN80, prop_postN30, prop_postE, prop_postG30, prop_postG80);

clear file_string ID session choice_criterion avgbet bets_10 bets_50 ...
    medianrt_bet10postloss medianrt_bet10postgain medianrt_bet50postloss medianrt_bet50postgain ...
    varrt_bet10postloss varrt_bet10postgain varrt_bet50postloss varrt_bet50postgain ...
    avgbet_postP10 avgbet_postP50 avgbet_postN10 avgbet_postN50 ...
    stdbet_postP10 stdbet_postP50 stdbet_postN10 stdbet_postN50 ...
    crit_postP10 crit_postP50 crit_postN10 crit_postN50 ...
    bets_postN80 bets_postN30 bets_postE bets_postG30 bets_postG80 ...
    avgbet_postN80 avgbet_postN30 avgbet_postE avgbet_postG30 avgbet_postG80 ...
    stdbet_postN80 stdbet_postN30 stdbet_postE stdbet_postG30 stdbet_postG80 ...
    prop prop_postP10 prop_postP50 prop_postN10 prop_postN50 ...
    prop_postN80 prop_postN30 prop_postE prop_postG30 prop_postG80

%%

writetable(beh_table, 'mf90_behav.csv');