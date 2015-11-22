function behdata=behav_fixnames(behdata)
% fix names in behdata

behdata.avgbet_po=behdata.avgbet_prevoutcome;
behdata.stdbet_po=behdata.stdbet_prevoutcome;

behdata.crit_po=behdata.crit_prevoutcome;

behdata.po_mat=behdata.prevoutcome_mat;

behdata.avgbet_po2=behdata.avgbet2_prevoutcome;
behdata.stdbet_po2=behdata.stdbet2_prevoutcome;

fields={'avgbet_prevoutcome','stdbet_prevoutcome','crit_prevoutcome', ...
    'prevoutcome_mat','avgbet2_prevoutcome','stdbet2_prevoutcome'};
behdata=rmfield(behdata,fields);

end