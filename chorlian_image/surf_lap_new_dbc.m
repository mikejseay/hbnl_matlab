function [fit,fit2,sig,snr,num]=surf_lap_new_dbc(x,y,xx,order);
%%% this program computes local quadratic fitting for random design and 
%%% for any output grid. If design and output grid are the same 2D
%%% redular grid, then loc_fit0.m is fast!!!
%%%  
%%% x: design matrix
%%% y: data
%%% xx: out_put grid
%%% order: polynomial order (2 or 3)
%%% USAGE
%%% 	[fit,fit2,sig,snr,num]=surf_lap(x,y,xx,order);
%%% fit: fitted surface
%%% fit2: estimated laplacian
%%% sig: estimated noise level (standard deviation)
%%% snr: estimated SNR
%%% num: number of neighboring points

[n,m]=size(x);
[nn,m]=size(xx);
fit=zeros(nn,1);
fit2=fit;

%%% u_b=10;
%%% l_b=10;
%%% n_min=15;
n_max=fix(n/2);
%%% for i=1:6,
%%% 	snr0(i)=l_b;
%%% 	l_b=l_b/2;
%%% 	nb(i)=round(n_min*(u_b/l_b)^(1/5));
%%% end;
%%% snr0=[100,snr0,0];
%%% ii=length(snr0);
%%% nb=[n_min,nb];

%%%%%%%%%%%%%%%%%%%%%
l_b=0.1;
nb_b=11;
for i=1:round(n/2)-11,
	nb(i)=nb_b+i;
	sig0(i)=l_b*(nb(i)/nb_b)^(9/2);
end;
sig0=[0,sig0,1000];
nb=[nb_b,nb];
iii=length(sig0);
%%%%%%%%%%%%%%%%%%%%%


%%%% estimate noise level sig
[ff,sig]=noise_est(x,y);
%%% ff_p=amse_error(ff);
%%% snr=ff_p/sig^2;
snr=max(0,(amse_error(y)-sig^2)/sig^2);

%%%% estimate bandwidth
%%% num=n_min;
%%% for i=1:ii-1,
%%% 	if (snr < snr0(i) & snr >= snr0(i+1)),
%%% 		num=min(nb(i),n_max);
%%% 	end;
%%% end;

%%%%%%%%%%%%%%%%%%%%%
num=nb_b;
for i=1:iii-1,
	if (sig < sig0(i+1) & sig >= sig0(i)),
		num=min(nb(i),n_max);
	end;
end;
%%%%%%%%%%%%%%%%%%%%%

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
% 	if (max(isnan(z)) == 1 | z(3) == 0),
% 		fit(i)=NaN;
% 		fit2(i)=NaN;
% 	else,
	
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
		
		[co,zzz,ind]=proj_tang(xz,z,10);
		
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
				%%%% w=(1-d1*d1);
				ssc=ssc+w*A*A';
				bc(:,j)=w*A;
			
			end;
		end;
	
		if rank(ssc) == mm,
			iss=inv(ssc)*bc;
		else,
	
			%%%% iss=inv(ssc+eps*EE)*bc;
			%%%% iss(1,:)=iss(1,:)/sum(iss(1,:));
			
			iss=inv(ssc)*bc;

		end;
		
		a=iss*yy;
	
		fit(i)=a(1);
		
		if mm > 3,
			fit2(i)=a(4)+a(5);
		end;
		
%	end;

end;
fit2 = -fit2;
