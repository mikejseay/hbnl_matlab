function out=meanx(in,except_dims)

%take a mean along all dimensions except those specified in except_dims
%ex: a is 4 x 9 x 36
% b = meanx(a,3) %returns a 1 x 36 vector of means along the other 2 dims

if nargin < 2
    except_dims = [];
end

n_dims=length(size(in));
if any(except_dims > n_dims)
    error('input matrix has fewer dimensions than specified')
end

mean_dims=setdiff(1:n_dims,except_dims);

out=in;
clear in
for d=mean_dims
    out=mean(out,d);
end
out=squeeze(out);

end