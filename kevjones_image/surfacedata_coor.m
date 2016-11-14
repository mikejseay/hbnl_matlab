function [surfdata] = surfacedata_coor(data, coor_new, n_points, surftype)
%
% surfacedata_coor [surfdata] =	surfacedata_coor(data, coor_new, n_points, surftype)
%
%			data is 2d data matrix (channels, time)
%			coor_new is x y and z channel positions
%			n_points is number of points for surfacing
%			surftype is surface type (0=surface potential [default], 1=surface laplacian)
%
%
if (nargin < 4)
	surftype = 0;
end

n_chans = size(coor_new, 1);
n_samples = size(data,2);

% load electrode coordinates

coor = coor_new;

nn = floor(n_points/2);
n=(2*nn+1)*(2*nn+1);
x_axis=(-nn:nn)/nn;
x=zeros(n,2);

for i=-nn:nn,
	for j=-nn:nn,
		x((i+nn)*(2*nn+1)+j+nn+1,:)=[i/nn,j/nn];
	end;
end;

[n,m]=size(x);
if m == 2, 
	x=[x,zeros(n,1)];
	for i=1:n,
		if norm(x(i,1:2)) > 1,
			x(i,3)=NaN;
		else,
			x(i,3)=sqrt(1-x(i,1:2)*x(i,1:2)');
		end;
	end;
end;

xx=x;

for i = 1:n_samples

	fprintf('sample : %d of %d \n', i, n_samples)
	y = data(:,i);
	yy=y(1:n_chans);
	
	[fit,fit2,sig,snr,num]=surf_lap_new_kevjones(coor(1:n_chans,:),yy,xx,2);	
	
	xysize = floor(sqrt(n));

	fitsurf = (reshape(fit,[xysize,xysize])');
	fitlap = (reshape(fit2,[xysize,xysize])');
	
	if (surftype == 1)
		n_points = size(fitlap,1);
		for j=1:n_points
			for k=1:n_points
				surfdata(i,j,k) = fitlap(j,k) * -1.0;
			end
		end
	else
		n_points = size(fitsurf,1);
		for j=1:n_points
			for k=1:n_points
				surfdata(i,j,k) = fitsurf(j,k);
			end
		end	
	end
	
end	
