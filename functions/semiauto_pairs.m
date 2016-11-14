% quick script to check the variance across each regional hypothesis' pair
% amounts

%% params

load('/export/home/mike/matlab/origin/coords/61chans_ns.mat')
type = 'regional';
% regions = {'frontal','central','parietal','occipital','lefttemporal','righttemporal'};
regions = {'frontal','parietal','occipital'};
min_dist = 5;
maxpairs2test = 12:46;
minpairs2test = 7;

%% generate initial pair set

if strcmpi(type,'combo') % combinatorial start point
    
    p = nchoosek(1:61,2);
    pi = 1:length(p);
    
elseif strcmpi(type,'regional') % regional start point

    [p,pi,pi_lbl] = determine_regional_pairs(chan_locs, 'both', regions, 2);

end

%% filter by distance and angle, remove tiny groups

[p2,pi2]=pair_filter(p,pi,chan_locs,'distance',[min_dist Inf]);
[p3,pi3]=pair_filter(p2,pi2,chan_locs,'angle','aplr');

n_hyps = max(pi3);
p4=[]; pi4=[]; pi_lbl2=cell(1,1); hd=0;
for hyp=1:n_hyps
    if sum(pi3==hyp) >= minpairs2test
        hd=hd+1;
        p4 = [p4; p3(pi3==hyp,:)];
        pi4 = [pi4; hd*ones(size(pi3(pi3==hyp)))];
        pi_lbl2{hd} = pi_lbl{hyp};
    end
end

%% test max pair amounts

n_hyps = max(pi4);
n_mpvals = length(maxpairs2test);
sp_dims = numSubplots(n_mpvals);
figure;
spd=0;
for pair_thresh = maxpairs2test
    spd=spd+1;
    subplot(sp_dims(1), sp_dims(2), spd);
    
    [p5,pi5] = sparsify_regional_hyps( p4, pi4, chan_locs, pair_thresh);
    
    [co,ce] = hist(pi5,n_hyps);
    pairvar(pair_thresh) = var(co);
    
    bar(ce,co);
    axis([0 n_hyps 0 maxpairs2test(end)]);
    title(['MPPH=',num2str(pair_thresh),', var=',num2str(pairvar(pair_thresh),3)])
    
end
figure;
subplot(211); plot(maxpairs2test,pairvar(maxpairs2test));
subplot(212); plot(maxpairs2test(2:end),diff(pairvar(maxpairs2test)));


%% pick a final amount and re-label for saving

[p5,pi5] = sparsify_regional_hyps( p4, pi4, chan_locs, 23);
figure; hist(pi5, 6);

%%

hyp_inds = pi5;
pairs = p5;
hyp_indlbls = pi_lbl2;

