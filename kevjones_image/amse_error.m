function ee=amse_error(x,y);

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
