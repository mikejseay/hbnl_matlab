function cmap=makecmap(v)
% makes a color map with white aligned with 0 and colorbrewer gradients
% toward positive (red) and negative limits (blue) indicated by 2-element
% vector v. if unsigned, uses a gradient from low (white) to high (green).

if v(1)<0 && v(2)>0 %(signed data)
    
    [~,zero_ind]=min(abs(linspace(v(1),v(2),256)));
    if 256-zero_ind > zero_ind
        cmap_pos=cbrewer('seq','Reds',256-zero_ind);
        cmap_neg=flipud(cbrewer('seq','Blues',256-zero_ind));
        cmap_neg(1:256-2*zero_ind,:)=[];
        cmap=[cmap_neg;cmap_pos];
    else
        cmap_pos=cbrewer('seq','Reds',zero_ind);
        cmap_neg=flipud(cbrewer('seq','Blues',zero_ind));
        cmap_pos(257-zero_ind:end,:)=[];
        cmap=[cmap_neg;cmap_pos];
    end
    
else
       
    cmap=cbrewer('seq','Greens',256);
    
end
        



end