function [estimates, model] = exp_fit4(data, dist, start_point, weight)
	% data ~ exp(a - b dist)
	if nargin < 3
		start_point = [0.8, 1, -2, 0.2];
	end
	if nargin < 4
		weight = ones(size(dist));
	end
	model = @exp_fit_model;
	estimates = fminsearch(model, start_point);
	
	function [errs, est] = exp_fit_model(params)
		a = real(params(1));
		b = real(params(2));
        c = real(params(3));
        d = real(params(4));
		est = a .^ (b * (dist + c) ) + d;
		err_vec = (data - est) .* weight;
		errs = sum(err_vec .^ 2);
	end
end
