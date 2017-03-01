% understand plotting phase as contourf

tf_scheme='cubicl'; %'linlhot' is also good

a=squeeze(angle(wavelet_evk_all(:,25,:,1,1)));
b=squeeze(real(wavelet_evk_all(:,25,:,1,1)));

%% contour of power-weighted phase in TF

cmap=pmkmp(256,tf_scheme);


figure
%imagesc(b')
contourf(fliplr(b)')
shading flat; colormap(cmap);

%% imagesc of phase in TF

cmap=pmkmp(256,tf_scheme);

figure
%imagesc(rad2deg(a)')
imagesc(rad2deg(a)')
caxis([-180 180]);
cmap=[cmap;flipud(cmap)];
colormap(cmap);

%% contour of cosine of phase in TF (eliminates contour edge effects)

cmap=pmkmp(256,tf_scheme);

figure
contourf(fliplr(cos(a))')
caxis([-1 1]);
shading flat; colormap(cmap);

%% contour of transformed phase in TF (eliminates contour edge effects)

c=squeeze(angle(wavelet_evk_all(:,25,:,1,1)));

cmap=pmkmp(256,tf_scheme);

c=pi-abs(2*pi-2*(c+pi));
figure
contourf(fliplr(c)')
shading flat; colormap(cmap);
%caxis([0 2*pi])

%% contour of transformed phase in TF (eliminates contour edge effects)
%in degrees

c=rad2deg(squeeze(angle(wavelet_evk_all(:,25,:,1,1)))+pi);

cmap=pmkmp(256,tf_scheme);

c=360-abs(360-c*2);
figure
contourf(fliplr(c)')
shading flat; colormap(cmap);
caxis([0 360])