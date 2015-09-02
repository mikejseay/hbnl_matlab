%  mobilize behavioral data into better matrices

ns=length(behdata);
n_conds=size(behdata(1).respRT,2);

meanRT=zeros(ns,n_conds);
medianRT=zeros(ns,n_conds);
varRT=zeros(ns,n_conds);
stdRT=zeros(ns,n_conds);

for s=1:ns
    
    meanRT(s,:)=behdata(s).respRT(1,:);
    medianRT(s,:)=behdata(s).respRT(2,:);
    varRT(s,:)=behdata(s).respRT(3,:);
    stdRT(s,:)=sqrt(behdata(s).respRT(3,:));
    
end