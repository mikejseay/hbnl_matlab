% testing various CSD applications

% full_load  % (load up something with a decent ERP)
load('/active_projects/matlab_common/61chans_CSD_GH.mat')

chan_erps = squeeze(erpdata(:,:,1,:));
coords = [chan_locs(:).X; chan_locs(:).Y; chan_locs(:).Z]';

%% original data
% is timepoints x channels

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-7 10]);
    set_print_size(18, 8);
end

%% CSD toolbox (spherical spline-based)

% lamda = 4
chan_erps_tb4 = zeros(size(chan_erps));

for s=1:size(chan_erps_tb4, 3)
    chan_erps_tb4(:, :, s) = CSD(squeeze(chan_erps(:,:,s))', csd_G, csd_H, 1.0e-4, 10)';
end

%%

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps_tb4(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps_tb4(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-.3 .3]);
    set_print_size(18, 8);
end

%% after mean

chan_erps_tb4_agmn = CSD(squeeze(mean(chan_erps,3))', csd_G, csd_H, 1.0e-4, 10)';

figure;
subplot(121); plot(chan_erps_tb4_agmn);
vline(248);
subplot(122); topoplot(chan_erps_tb4_agmn(248,:),chan_locs);
set_print_size(18, 8);

%%

figure;
subplot(121); plot(chan_erps_tb4_gmn);
vline(248);
subplot(122); topoplot(chan_erps_tb4(248,:),chan_locs);
set_print_size(18, 8);

%%
% lamda = 5
chan_erps_tb5 = zeros(size(chan_erps));

for s=1:size(chan_erps_tb5, 3)
    chan_erps_tb5(:, :, s) = CSD(squeeze(chan_erps(:,:,s))', csd_G, csd_H, 1.0e-5, 10)';
end

%%

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps_tb5(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps_tb5(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-.3 .3]);
    set_print_size(18, 8);
end

%%
% lamda = 5.5
chan_erps_tb55 = zeros(size(chan_erps));

for s=1:size(chan_erps_tb55, 3)
    chan_erps_tb55(:, :, s) = CSD(squeeze(chan_erps(:,:,s))', csd_G, csd_H, 5.0e-6, 10)';
end

%%

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps_tb55(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps_tb55(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-.3 .3]);
    set_print_size(18, 8);
end


%%
% lamda = 6
chan_erps_tb6 = zeros(size(chan_erps));

for s=1:size(chan_erps_tb6, 3)
    chan_erps_tb6(:, :, s) = CSD(squeeze(chan_erps(:,:,s))', csd_G, csd_H, 1.0e-6, 10)';
end

%%

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps_tb6(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps_tb6(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-.3 .3]);
    set_print_size(18, 8);
end


%% David Chorlian's rewrite of Kongming Wang's local polynomical

chan_erps_dbc = zeros(size(chan_erps));

for s=1:size(chan_erps_tb6, 3)
    chan_erps_dbc(:, :, s) = surf_laplacian_dbc(squeeze(chan_erps(:,:,s)), coords, 0, 11);
end

%%

for g = pp.chosen_g(pp.plotn_g)
    figure;
    subplot(121); plot(squeeze(mean(chan_erps_dbc(:, :, s_inds_g(:,g)),3)));
    vline(248);
    subplot(122); topoplot(squeeze(mean(chan_erps_dbc(248, :, s_inds_g(:,g)),3)), chan_locs, ...
        'maplimits', [-30 30]);
    set_print_size(18, 8);
end

%% Kevin Jones' MATLAB-based implementation (based on Kongming Wang's local polynomial)

% slow, so we only do 1 point
chan_erps_kevmat = surface_channel_data(chan_erps(248,:), coords, 150, 1);

%%

figure; contourf(squeeze(chan_erps_kevmat(1, :, :)));
set(gca, 'xdir', 'reverse');

%% Kevin Jone's 2nd MATLAB-based implementation (also slow)

chan_erps_kevmat2 = surfacedata_coor(chan_erps(248,:)', coords, 150, 1);

%%

figure; contourf(squeeze(chan_erps_kevmat2(1, :, :)));
set(gca, 'xdir', 'reverse');

%% Kevin Jones' C-based implementation
% in 'surf_hdf1_laplace' which can only be run on v490

path = '/export/home/mike/avg_hp003_t75/vp3_5_a1_a0001410_csd.h1';
h1_struct = read_hdf1_dataV7(path);

chan_vec = [1:31,33:62];
sen = 204.8./6.1035156;
chan_erps_kevc = squeeze(h1_struct.data_struct.hdf1_avg_data(1,chan_vec,:))'*sen;

figure;
subplot(121); plot(chan_erps_kevc);
vline(248);
subplot(122); topoplot(chan_erps_kevc(168,:),chan_locs);
set_print_size(18, 8);

%% EEGLAB's function (seems very off)

%{
chan_topo_lab = del2map(chan_erps(248,:)', '/active_projects/matlab_common/61chans_ns.ced')';
figure; topoplot(chan_topo_lab ,chan_locs);
%}
