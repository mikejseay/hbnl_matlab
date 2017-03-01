% assess differences between non-cleaned vs. cleaned data

chan=7;
scale=15;
cond=1;

ero_data=Y.wavelet_tot;
ero_data_old=Y_old.wavelet_tot;

itc_data=abs(Y.wavelet_evk)./Y.wavelet_tot;
itc_data_old=abs(Y_old.wavelet_evk)./Y_old.wavelet_tot;

ero_image_data=squeeze(ero_data(:,chan,:,cond));
ero_image_data_old=squeeze(ero_data_old(:,chan,:,cond));

itc_image_data=squeeze(itc_data(:,chan,:,cond));
itc_image_data_old=squeeze(itc_data_old(:,chan,:,cond));

ero_lims=[5 60];
itc_lims=[0 0.8];

%%

figure;
subplot(2,1,1);
contourf(fliplr(ero_image_data)'); shading flat;
caxis(ero_lims)
subplot(2,1,2);
contourf(fliplr(ero_image_data_old)'); shading flat;
caxis(ero_lims)
colormap(pmkmp(256,'cubicl'));
title('ERO')

%%

figure;
subplot(2,1,1);
contourf(fliplr(itc_image_data)'); shading flat;
caxis(itc_lims)
subplot(2,1,2);
contourf(fliplr(itc_image_data_old)'); shading flat;
caxis(itc_lims)
colormap(pmkmp(256,'cubicl'));
title('ITC')

%%

fprintf('Cleaning yielded an average of %1.1f more trials per condition',...
    round(mean(Y.n_trials-Y_old.n_trials)))