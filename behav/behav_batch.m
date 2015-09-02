%quick script to check the idea of bet amount being based on previous
%outcome

init_type=90;

n_files=length(h1_list);
%avgprev=zeros(n_files,4);
%stdprev=zeros(n_files,4);
%crit=zeros(n_files,1);
%crit_prevoutcome=zeros(n_files,4);
%avgbet=zeros(n_files,1);
%n_bets=zeros(n_files,2);
%combout=zeros(n_files,5);
%betaftertrend=zeros(n_files,5);
%avgprev2=zeros(n_files,5);
%stdprev2=zeros(n_files,5);

for file=1:n_files
    tic
    h1_struct=read_hdf1_behdata(h1_list{file});
    etable=h1_getbehav(h1_struct,init_type);
    [avgprev2(file,:),stdprev2(file,:)]=behav_sog(etable);
    toc
end