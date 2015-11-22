% speed / space testing for data formats

time = [];

%% mat

tic;
save(output_file, '-struct', 'output', '-v6');
time.matv6 = toc;

%% h5

tic;
h5create(output_file,filename(1,:),size(output.wave_tot),'Datatype','double');
h5write(output_file,filename(1,:),output.wave_tot);
time.h5 = toc;

%% json