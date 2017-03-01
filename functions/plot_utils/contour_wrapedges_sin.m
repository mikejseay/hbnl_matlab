function [data_out] = contour_wrapedges_sin(data_in)

data_out=-pi+abs(2*pi-2*(data_in+pi));

end