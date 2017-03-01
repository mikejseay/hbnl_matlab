% testing various CSD applications

% full_load  % (load up something with a decent ERP)
load('/active_projects/matlab_common/61chans_CSD_GH.mat')

chan_erps = squeeze(erpdata(:,:,1,3));
coords = [chan_locs(:).X; chan_locs(:).Y; chan_locs(:).Z]';

%% original data
% is timepoints x channels

figure;
subplot(121); plot(chan_erps);
vline(248);
subplot(122); topoplot(chan_erps(248,:),chan_locs);
set_print_size(18, 8);

%% CSD toolbox (spherical spline-based)

% lamda = 4

chan_erps_tb4 = CSD(chan_erps', csd_G, csd_H, 1.0e-4, 10)';

figure;
subplot(121); plot(chan_erps_tb4);
vline(248);
subplot(122); topoplot(chan_erps_tb4(248,:),chan_locs);
set_print_size(18, 8);

%%

% lambda = 5 (default)

chan_erps_tb = CSD(chan_erps', csd_G, csd_H, 1.0e-5, 10)';

figure;
subplot(121); plot(chan_erps_tb);
vline(248);
subplot(122); topoplot(chan_erps_tb(248,:),chan_locs);
set_print_size(18, 8);

%%

% lambda = 5.5

chan_erps_tb = CSD(chan_erps', csd_G, csd_H, 5.0e-6, 10)';

figure;
subplot(121); plot(chan_erps_tb);
vline(248);
subplot(122); topoplot(chan_erps_tb(248,:),chan_locs);
set_print_size(18, 8);

%%

% lambda = 6

chan_erps_tb6 = CSD(chan_erps', csd_G, csd_H, 1.0e-6, 10)';

figure;
subplot(121); plot(chan_erps_tb6);
vline(248);
subplot(122); topoplot(chan_erps_tb6(248,:),chan_locs);
set_print_size(18, 8);

%% David Chorlian's rewrite of Kongming Wang's local polynomical

chan_erps_dbc = surf_laplacian_dbc(chan_erps, coords, 0, 11);

figure;
subplot(121); plot(chan_erps_dbc);
vline(248);
subplot(122); topoplot(chan_erps_dbc(248,:),chan_locs);
set_print_size(18, 8);

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
