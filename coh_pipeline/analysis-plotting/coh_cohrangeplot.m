% wonderful idea: i wanna be like d

%% calculate the 'chan_cohlims' variable

chan_cohlims=meanx(cohlim_mat,[1 4])';

%% each pair's "coherence range" as a line from min (left) to max (right)
figure;
for pair=1:size(chan_cohlims,1)
plot(chan_cohlims(pair,:),'Color',scl.p_color(pair,:)); hold on;
end
axis([0 3 0 1])

% color line by distance of pair

%% each pairs "coherence range" as a vertical value for each pair
chan_cohranges=diff(chan_cohlims,1,2);
figure; plot(chan_cohranges);

%% scattering pair distance (x axis) against minimum / maximum coherence
p_dists=pair_distance(opt.coherence_pairs,chan_locs);
figure; scatter(p_dists,chan_cohlims(:,1),'b'); hold on;
scatter(p_dists,chan_cohlims(:,2),'r');

%% plotting each range by distance
p_dists=pair_distance(opt.coherence_pairs,chan_locs);
figure;
for pair=1:size(chan_cohlims,1)
    line([p_dists(pair) p_dists(pair)],chan_cohlims(pair,:),'Marker','+','Color',scl.p_color(pair,:)); hold on;
end
ylabel('Coherence Min/Max')
xlabel('Distance (cm)')

%% scattering pair distance against coherence range
figure; scatter(p_dists,chan_cohranges);
xlabel('Distance (cm)')
ylabel('Coherence Range')