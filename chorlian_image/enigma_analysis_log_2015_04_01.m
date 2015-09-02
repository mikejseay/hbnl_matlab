clear all
file = '/processed_data/eeg/enigma.d/ec8A0006001.rd.bin';
load eeg_monopolar_struct.mat
X = load('eeg_bipolar_struct.mat');
X = load('eeg_mono_structB.mat');
X = load('eeg_mono_structC.mat');
[data, count] = read_bin(file, 256 *256, 19);
file = '/processed_data/eeg/enigma.d/data/ec8A0006001.rd.bin';
[data, count] = read_bin(file, 256 *256, 19);
figure; plot(data);
figure; plot(data(:, 1:3))
figure; plot(data(:, 1))
hold on; plot(data(:, 2) + 500)
hold on; plot(data(:, 3) - 500)
figure; for m = 1:10
plot(data(:, m) - (m -1) * 500);
hold on
end; grid on
figure; plot(data(:,1) + 500); hold on; plot(data(:, 15)); grid on
data_ec8A0006001 = data;
file = '/processed_data/eeg/enigma.d/data/ec_8_a1_50501001.cnt.h5.rd.bin';
[data, count] = read_bin(file, 256 *256, 20);
file = '/processed_data/eeg/enigma.d/data/ec8A0006001.rd.bin';
data_ec8A0006001 =  read_bin(file, 256 *256, 20);
figure; plot(data(:,1) + 500); hold on; plot(data(:, 15)); grid on
data_50501001 = data;
file = '/processed_data/eeg/enigma.d/data/ec8a0071001.rd.bin';
[data, count] = read_bin(file, 256 *256, 20);
figure; plot(data(:,1) + 500); hold on; plot(data(:, 15)); grid on
std(data(:,1))
std(data)
std(data_ec8A0006001)
sum(data(:,1) == 0)
sum(data_50501001(:,1) == 0)
figure; plot(data_50501001(1:1024,1)); grid on
figure; plot(data_50501001(1:10240,1)); grid on
figure; plot(data_50501001(1024:4096,1)); grid on
figure; plot(data_50501001(2048:4096,1)); grid on
open('enigma_pre_analysis.m');
file = '/processed_data/eeg/enigma.d/data/ec_8_a1_50501001.cnt.h5.rd.bin';
[Y_50501001, T_50501001] = enigma_pre_analysis(file);
figure; plot(double(T_50501001), 'x');
open('enigma_analysis.m');
dbstop at 1 in enigma_analysis
enigma_analysis(file);
dbquit
dbstop at 1 in enigma_analysis
enigma_analysis(file);
dbstop if error
dbstop at 1 in enigma_analysis
enigma_analysis(file);
dbquit
enigma_analysis(file);
size(Y(trial_start(trial):trial_finish(trial + 1), :))
size(hanning_window)
dbquit
dbstop at 1 in enigma_analysis
enigma_analysis(file);
figure; plot(P)
figure; plot(PW)
axis([0, 64, 0, 6]); grid on
figure; plot(fs, PW); grid on
axis([0, 40, 0, 6]); grid on
m = 1; f_idx = fs >= cFreq(m, 1) & fs <=  cFreq(m, 2);
PF = mean(log(PW(fs, :)));
PF = mean(log(PW(f_idx, :)));
size(PF)
PF
log(mean(PW((f_idx, :)))
log(mean(PW(f_idx, :)))
mean(PW(f_idx, :))
X = log(PW(f_idx, :));
XX = PW(f_idx, :);
f_idx
dbquit
file = '/processed_data/eeg/ns_32_eeg/eec_3_a1_30125020_32.cnt.h1';
h1_struct = read_hdf1_dataV7(data_file);
h1_struct = read_hdf1_dataV7(file);
X = h1_struct.data_struct.hdf1_cnt_data;
data = X';
clear X
h1_chan_idx = [16, 31, 30];
X = data(:, h1_chan_idx);
figure; plot(X(:, 1));
XX = detrend(X);
figure; plot(XX(:, 1));
dataF = filter_data(XX, .2, 40, [], 256);
figure; plot(dataF(:, 1))
axis([0, 1024, -Inf, Inf]);
hanning_window = repmat(hanning(512), 1, 3);
Fh =  hanning_window .* dataF(1:512, :);
figure; plot(Fh);
figure; plot(Fh(:, 1))
[Xd, Xt] = enigma_artf_check(X);
[XFd, XFt] = enigma_artf_check(dataF);
sum(XFt ~= Xt)
fidx = XFt ~= X;
fidx = XFt ~= Xt;
XXX = [Xt, XFt];
Z = XXX(fidx);
Z = XXX(fidx, :);
dbclear all
open('enigma_analysis.m');
dbstop at 1 in enigma_analysis
enigma_analysis(XFd, XFt);
dbquit
enigma_analysis(XFd, XFt);
sum(trial_idx)
dbquit
sum(Xt)
file = '/processed_data/eeg/filtered_data/coga/iowa/ec_8_a1_30125020_32.cnt.h1.rd';
dbstop at 1 in  enigma_analysis_mc
enigma_analysis_mc(file);
dbquit
file = '/processed_data/eeg/filtered_data/coga/washu/ec_8_a1_50019058_32.cnt.h1.rd';
enigma_analysis_mc(file);
dbstop at 1 in  enigma_analysis_mc
enigma_analysis_mc(file);
dbclear all
file = '/processed_data/eeg/filtered_data/coga/washu/ec_8_a1_50071010_32.cnt.h1.rd';
output = enigma_analysis_mc(file);
file = '/processed_data/eeg/filtered_data/coga/washu/ec_8_a1_50189023_32.cnt.h1.rd';
output = enigma_analysis_mc(file);
