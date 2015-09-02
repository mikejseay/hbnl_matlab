% understand RMANOVA in matlab

load fisheriris

t = table(species,meas(:,1),meas(:,2),meas(:,3),meas(:,4),...
'VariableNames',{'species','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Measurements'});

rm = fitrm(t,'meas1-meas4~species','WithinDesign',Meas);
rm_result= ranova(rm);