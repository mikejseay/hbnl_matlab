function [co,y,ind]=proj_tang(x,x0,h);
%%% x: points on unit sphere
%%% x0: a point on unit sphere
%%% h: bandwidth 
%%% co: coordinates on the local surface
%%% y: projection of x on the tangent plan at x0
%%% ind: index of data points with the local coordinates

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
	
	
	%%%%%%%%%%%%%%%%% keep the surface distance
	%%%%% if norm(y(i,:)-x0) > 0;
	%%%%% 	dd=(1-lambda)/norm(y(i,:)-x0);
	%%%%% else,
	%%%%% 	dd=0;
	%%%%% end;
	%%%%% co(i,1)=dd*(y(i,:)-x0)*coor1';
	%%%%% co(i,2)=dd*(y(i,:)-x0)*coor2';
	%%%%%%%%%%%%%%%%% 
	
	%%%%%%%%%%%%%%%%% distance decreased
	co(i,1)=(y(i,:)-x0)*coor1';
	co(i,2)=(y(i,:)-x0)*coor2';
	%%%%%%%%%%%%%%%%% 

end;

