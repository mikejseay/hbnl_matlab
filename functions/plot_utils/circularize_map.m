function cmap_out=circularize_map(cmap_in)
% dumb simple way of circularizing any color map - will not work for
% multi-chromatic colors but will work for linear monochromatic maps

%all channels simultaneously

%record the incoming minimum and maximum values
inminval=min(min(cmap_in));
inmaxval=max(max(cmap_in));

% normalize to [0 1]
cmap_norm=norm2limits_mat(cmap_in,[0 1]);

%transform to cosine
cmap_cos=cos(cmap_norm*2*pi);

%scale back to the in min and in max
cmap_out=norm2limits_mat(cmap_cos,[inminval inmaxval]);

%normalize back to [0 1] for each channel separately

%{
for channel=1:3
    
    minval=min(cmap_cos(:,channel)); %could be as low as -1

    cmap_pos(:,channel) = ( cmap_cos(:,channel) + abs(minval) ); % all positive vals now
   
    maxval=max( cmap_pos(:,channel) ); %could be as high as 2 i guess

    cmap_out(:,channel) = cmap_pos(:,channel) / maxval; %now 0 to 1
    
end
%}

end