function [surfdata] = surface_channel_data(data, coor, n_points, surftype)
%
% [surfdata] =	surface_channel_data(data, coor, n_points, surftype)
%
%			data is 2d data matrix (time by channels)
%			coor is x y and z channel positions
%			n_points is number of points for surfacing
%			surftype is surface type (0=surface potential [default], 1=surface laplacian)
%
%
	if (nargin < 4)
		surftype = 0;
	end

	n_chans = size(coor,1);
	n_samples = size(data,1);

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

	% setup laplacian design matrix

	y = data(1,:);
	yy=y(1:n_chans)';
	[lap_design_matrix,num]=setup_surface_laplacian(coor(1:n_chans,:),yy,xx,2);

	% apply laplacian surfacing

	for i = 1:n_samples

		y = data(i,:);
		yy=y(1:n_chans)';

		%[fit,fit2,sig,snr,num]=do_surface_laplacian(coor(1:n_chans,:),yy,xx,2);	

		[fit,fit2]=apply_surface_laplacian(coor(1:n_chans,:),yy,xx,lap_design_matrix,num);

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

function [lap_design_matrix,num]=setup_surface_laplacian(x,y,xx,order);
%
% function [lap_design_matrix,num]=setup_surface_laplacian(x,y,xx,order);
%
%%% this program computes local quadratic fitting for random design and 
%%% for any output grid. If design and output grid are the same 2D
%%% redular grid, then loc_fit0.m is fast
%%%  
%%% x: design matrix
%%% y: data
%%% xx: out_put grid
%%% order: polynomial order (2 or 3)
%%% 
%%% lap_design_matrix: laplacian design matrix
%%% num
%%%

	[n,m]=size(x);
	[nn,m]=size(xx);
	fit=zeros(nn,1);
	fit2=fit;

	n_max=fix(n/2);

	nb=[];
	sig0=[];
	l_b=0.1;
	nb_b=11;
	for i=1:round(n/2)-11,
		nb(i)=nb_b+i;
		sig0(i)=l_b*(nb(i)/nb_b)^(9/2);
	end;
	sig0=[0,sig0,1000];
	nb=[nb_b,nb];
	iii=length(sig0);

	%%%% estimate noise level sig

	[ff,sig]=noise_estimate(x,y);
	snr=max(0,(asymptotic_mean_square_error(y)-sig^2)/sig^2);

	num=nb_b;
	
	for i=1:iii-1,
		if (sig < sig0(i+1) & sig >= sig0(i)),
			num=min(nb(i),n_max);
		end;
	end;

	if order == 2,
		mm=6;   %%% local quadratic fitting
	elseif order == 3,
		mm=8;
	else,
		mm=6;
	end;

	EE=eye(mm,mm);
	eps=0.0001;
	dd=zeros(n,1);

	for i=1:nn,

		z=xx(i,:);
		if (max(isnan(z)) == 1 | z(3) == 0),
			fit(i)=NaN;
			fit2(i)=NaN;
		else,

			for j=1:n,
				dd(j)=norm(x(j,:)-z);
			end;
			[dd,ind]=sort(dd);

			nm=num;
			for kk=num+1:n,
				if dd(kk) == dd(num),
					nm=kk;
				end;
			end;
			xz=x(ind(1:nm),:);
			yy=y(ind(1:nm));
			dddd=(dd(nm+1)-dd(nm))/2;

			[co,zzz,ind]=project_tangent(xz,z,10);

			for j=1:nm,
				dd(j)=norm(co(j,:));
			end;

			ssc=zeros(mm,mm);
			bc=zeros(mm,nm);
			h=max(dd(1:nm))+dddd;

			for j=1:nm,

				if dd(j) > h,
					bc(:,j)=zeros(mm,1);
				else,
					if mm == 6,
						A=[1;co(j,:)';(co(j,:).^2/2)';co(j,1)*co(j,2)]; 
					elseif mm == 8,
						A=[1;co(j,:)';(co(j,:).^2/2)';co(j,1)*co(j,2);(co(j,:).^3/6)'];
					else,
						A=[1;co(j,:)'];
					end;

					d1=dd(j)/h;
					w=(1-d1*d1)^2;  % kernel
					ssc=ssc+w*A*A';
					bc(:,j)=w*A;

				end;
			end;

			if rank(ssc) == mm,
				iss=inv(ssc)*bc;
			else,			
				iss=inv(ssc)*bc;
			end;

			lap_design_matrix{i} = iss;

		end;

	end;

function [fit,fit2]=apply_surface_laplacian(x,y,xx,lap_design_matrix,num);
%
% function [fit,fit2]=apply_surface_laplacian(x,y,xx,lap_design_matrix,num);
%
%%%  
%%% x: design matrix
%%% y: data
%%% xx: out_put grid
%%% lap_design_matrix: pre-calculated design matrices
%%% num
%%%
%%% fit: fitted surface
%%% fit2: estimated laplacian
%%% 

	[n,m]=size(x);
	[nn,m]=size(xx);
	fit=zeros(nn,1);
	fit2=fit;
	n_max=fix(n/2);
	
	for i=1:nn,

		z=xx(i,:);

		if (max(isnan(z)) == 1 | z(3) == 0),
			fit(i)=NaN;
			fit2(i)=NaN;
		else,
			
			for j=1:n,
				dd(j)=norm(x(j,:)-z);
			end;
			[dd,ind]=sort(dd);

			nm=num;
			for kk=num+1:n,
				if dd(kk) == dd(num),
					nm=kk;
				end;
			end;

			yy=y(ind(1:nm));
			
			iss = lap_design_matrix{i};

			a=iss*yy;

			fit(i)=a(1);
			
			fit2(i)=a(4)+a(5);

		end;

	end;

function [fit,fit2,sig,snr,num]=do_surface_laplacian(x,y,xx,order);
%
% function [fit,fit2,sig,snr,num]=do_surface_laplacian(x,y,xx,order);
%
%%% this program computes local quadratic fitting for random design and 
%%% for any output grid. If design and output grid are the same 2D
%%% redular grid, then loc_fit0.m is fast!!!
%%%  
%%% x: design matrix
%%% y: data
%%% xx: out_put grid
%%% order: polynomial order (2 or 3)
%%% USAGE
%%% 	[fit,fit2,sig,snr,num]=surface_laplacian(x,y,xx,order);
%%% fit: fitted surface
%%% fit2: estimated laplacian
%%% sig: estimated noise level (standard deviation)
%%% snr: estimated SNR
%%% num: number of neighboring points
%%% 

	[n,m]=size(x);
	[nn,m]=size(xx);
	fit=zeros(nn,1);
	fit2=fit;

	n_max=fix(n/2);

	nb=[];
	sig0=[];
	l_b=0.1;
	nb_b=11;
	for i=1:round(n/2)-11,
		nb(i)=nb_b+i;
		sig0(i)=l_b*(nb(i)/nb_b)^(9/2);
	end;
	sig0=[0,sig0,1000];
	nb=[nb_b,nb];
	iii=length(sig0);

	%%%% estimate noise level sig

	[ff,sig]=noise_estimate(x,y);
	snr=max(0,(asymptotic_mean_square_error(y)-sig^2)/sig^2);

	num=nb_b;
	for i=1:iii-1,
		if (sig < sig0(i+1) & sig >= sig0(i)),
			num=min(nb(i),n_max);
		end;
	end;

	if order == 2,
		mm=6;   %%% local quadratic fitting
	elseif order == 3,
		mm=8;
	else,
		mm=3;
	end;

	EE=eye(mm,mm);
	eps=0.0001;
	dd=zeros(n,1);

	for i=1:nn,

		z=xx(i,:);
		if (max(isnan(z)) == 1 | z(3) == 0),
			fit(i)=NaN;
			fit2(i)=NaN;
		else,

			for j=1:n,
				dd(j)=norm(x(j,:)-z);
			end;
			[dd,ind]=sort(dd);

			nm=num;
			for kk=num+1:n,
				if dd(kk) == dd(num),
					nm=kk;
				end;
			end;
			xz=x(ind(1:nm),:);
			yy=y(ind(1:nm));
			dddd=(dd(nm+1)-dd(nm))/2;

			[co,zzz,ind]=project_tangent(xz,z,10);

			for j=1:nm,
				dd(j)=norm(co(j,:));
			end;

			ssc=zeros(mm,mm);
			bc=zeros(mm,nm);
			h=max(dd(1:nm))+dddd;

			for j=1:nm,

				if dd(j) > h,
					bc(:,j)=zeros(mm,1);
				else,
					if mm == 6,
						A=[1;co(j,:)';(co(j,:).^2/2)';co(j,1)*co(j,2)]; 
					elseif mm == 8,
						A=[1;co(j,:)';(co(j,:).^2/2)';co(j,1)*co(j,2);(co(j,:).^3/6)'];
					else,
						A=[1;co(j,:)'];
					end;

					d1=dd(j)/h;
					w=(1-d1*d1)^2;  % kernel
					ssc=ssc+w*A*A';
					bc(:,j)=w*A;

				end;
			end;

			if rank(ssc) == mm,
				iss=inv(ssc)*bc;
			else,			
				iss=inv(ssc)*bc;
			end;

			a=iss*yy;

			fit(i)=a(1);

			if mm > 3,
				fit2(i)=a(4)+a(5);
			end;

		end;

	end;

function ee=asymptotic_mean_square_error(x,y);

	[n,m]=size(x);

	if nargin < 2,
		y=zeros(n,m);
	end;

	z=x-y;

	k=0;
	for i=1:n,
		for j=1:m,
			if isnan(z(i,j)) == 1;
				z(i,j)=0;
				k=k+1;
			end;
		end;
	end;

	ee=sum(sum(z.*z))/(n*m-k);


function [ff,sig]=noise_estimate(x0,y0); 

	[n,m]=size(x0);

	yn=1-isnan(y0);

	x=[];
	y=[];

	for i=1:n,
		if yn(i) > 0,
			x=[x;x0(i,:)];
			y=[y;y0(i)];
		end;
	end;

	[n,m]=size(x);

	dd=zeros(n,1);
	ff=dd;

	for i=1:n,
		for j=1:n,
			dd(j)=norm(x(j,:)'-x(i,:)');
		end;
		[dd,ind]=sort(dd);

		nm=3;
		for j=5:n,
			if dd(j) == dd(4), nm=nm+1; end;
		end;

		xx=x(ind(2:nm+1),:);
		yy=y(ind(2:nm+1));

		[co,u,v]=project_tangent(xx,x(i,:),10);

		ff(i) = NaN;
		AA=[ones(nm,1),co];

		[uaa,saa,vaa]=svd(AA);

		if saa(3,3) < 0.0001,
			AINV=inv(AA'*AA + 0.0005*eye(3));
		else,
			AINV=inv(AA'*AA);
		end;

		a=AINV*(AA'*yy);

		ff(i)=a(1);
	end;

	sig=((y-ff)'*(y-ff)/(n-1)+median(abs(y-ff))^2)/2;

	sig=sqrt(sig);

function [co,y,ind]=project_tangent(x,x0,h);

%%% x: points on unit sphere
%%% x0: a point on unit sphere
%%% h: bandwidth 
%%% co: coordinates on the local surface
%%% y: projection of x on the tangent plan at x0
%%% ind: index of data points with the local coordinates
%%% 

	[n,m]=size(x);

	ind=[];

	k=0;
	for i=1:n,
		d=(x(i,:)-x0)*(x(i,:)-x0)';
		if d <= h*h, 
			k=k+1;
			ind=[ind,i];
		end;
	end;

	%%%% coordinate system on the local surface
	
	if x0(3) < 1,
		coor1=[-x0(1)*x0(3),-x0(2)*x0(3),1-x0(3)*x0(3)]/sqrt(1-x0(3)*x0(3));
		coor2=[x0(2),-x0(1),0]/sqrt(x0(1)*x0(1)+x0(2)*x0(2));
	else,
		coor1=[1,0,0];
		coor2=[0,1,0];
	end;

	%%%% projection of x on the tangent plan at x0 
	
	y=zeros(k,m);
	co=zeros(k,2);

	for i=1:k,

		lambda=1-x(ind(i),:)*x0';

		y(i,:)=x(ind(i),:)+lambda*x0;

		%%%% distance decreased
		co(i,1)=(y(i,:)-x0)*coor1';
		co(i,2)=(y(i,:)-x0)*coor2';

	end;
