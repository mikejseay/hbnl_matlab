function [ff,sig]=noise_est(x0,y0);

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

	[co,u,v]=proj_tang(xx,x(i,:),10);
	
	%%% a=inv([ones(nm,1),co])*yy;
	
	%%% if isempty(co) == 1,
		ff(i) = NaN;
	%%% else,
		AA=[ones(nm,1),co];
		
		[uaa,saa,vaa]=svd(AA);
		
		if saa(3,3) < 0.0001,
			AINV=inv(AA'*AA + 0.0005*eye(3));
		else,
			AINV=inv(AA'*AA);
		end;
		
		a=AINV*(AA'*yy);
	
		ff(i)=a(1);
	%%% end;
end;

%%% sig=(y-ff)'*(y-ff)/n;
%%% sig=median(abs(y-ff))^2;
sig=((y-ff)'*(y-ff)/(n-1)+median(abs(y-ff))^2)/2;

sig=sqrt(sig);
