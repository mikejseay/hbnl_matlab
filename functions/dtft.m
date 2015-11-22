function fourier=dtft(in,n_pts)

n = length(in);

if nargin<2
    n_pts=n;
end

if n_pts>n
    add_pts=n_pts - n;
    if mod(add_pts,2)==0
        in=[zeros(1,add_pts/2),in,zeros(1,add_pts/2)];
    else
        odd_pts=floor(add_pts/2);
        in=[zeros(1,odd_pts),in,zeros(1,odd_pts+1)]; %not sure
    end
end

fourier=zeros(size(in));
time=(0:n_pts-1)/n_pts;
for fi=1:n_pts
    sine=exp(-1i*2*pi*(fi-1).*time);
    fourier(fi) = sum(sine.*in);
end

end