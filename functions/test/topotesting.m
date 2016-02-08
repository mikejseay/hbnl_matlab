figure;

dim1=3;
dim2=2;

for sp=1:dim1*dim2
    testdata = randn([61 1]);
    sp_temp = subplot(dim1,dim2,sp);
    h=topoplot(testdata, chan_locs, 'electrodes', 'off', 'colormap', cmap, 'style', 'fill', 'numcontour', 7);
    if sp < 3
        colormap(sp_temp,cmap);
    else
        colormap(sp_temp,cmap_diff);
    end
    set(h,'EdgeColor','None');
end

% there is a threshold!!! it worked for 6 but not for 8!!