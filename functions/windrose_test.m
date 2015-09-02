% understand wind-rose

d=0:10:350;
D=[];
V=[];
for i=1:length(d)
n=d(i)/10;
D=[D ones(1,n)*d(i)];
V=[V 1:n];
end

figure
wind_rose(D,V)
%
figure
wind_rose(D,V,'iflip',1)
%
figure
wind_rose(D,V,'ci',[1 2 7],'dtype','meteo')

%%

angletest=squeeze(mean(angle(wavelet_evk_all(129:140,7,10,1,s_inds)),1));
%angletest_i=squeeze(mean(abs(wavelet_evk_all(129:140,7,10,1,s_inds)),1));
angletest_i=squeeze(mean(itc_results_all(129:140,7,10,1,s_inds),1));
angletest_r=rad2deg(angletest);
figure;
[handles,data]=wind_rose(angletest_r,angletest_i,'n',36,'di',[0:0.1:0.6],...
    'ci',3,'labtitle','title','lablegend','coherence','cmap',cmap,...
    'quad',3);
%figure; wind_rose(angletest_r,angletest_i);