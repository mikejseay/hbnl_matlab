function recon=idtft(in)

n_pts = length(in);

recon=zeros(size(in));
time=(0:n_pts-1)/n_pts;

in=in/n_pts; %"normalize" the fourier coefficients

for fi=1:n_pts
    sine=in(fi)*exp(1i*2*pi*(fi-1).*time);
    recon = recon + real(sine);
end

recon=circshift(recon,round(n_pts/2),2);

end