function out = norm2limits ( in, lims )
% normalize a single number to specified min/max limits

if lims(1) < 0 && lims(2) > 0 % i.e. polar data
    if in < 0
        out = abs( in / lims(1) );
    elseif in > 0
        out = abs( in / lims(2) );
    end
else
    out = ( in - lims(1) ) / range(lims); %i.e. non-polar data
end

end