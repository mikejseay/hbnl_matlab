function out = norm2limits_mat(in, lims)
% normalize a 2d matrix to provided limits

minval=min(min(in));

out1=in-minval; %lowest val is now 0

maxval=max(max(out1));

out2=out1./maxval; %highest val is now 1

out3 = out2 .* range(lims);

out = out3 + lims(1);

end