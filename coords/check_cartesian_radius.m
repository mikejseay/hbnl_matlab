radius=zeros(size(data,1),1);
for chan=1:size(data,1)
    radius(chan)=sqrt(data(chan,1)^2+data(chan,2)^2+data(chan,3)^2);
end