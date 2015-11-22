function [estimates, model] = exp_fit3(data, dist, start_point, weight)
	% data ~ exp(a + b dist) + c
	if nargin < 3 || isempty(start_point)
		b = regress(data, [ones(length(data), 1), exp(-dist)]);
		if b(2) < .001
			b(2) = .001;
		end
		start_point = [log(b(2)), -1, min(data)];
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
		est = exp(a + b * dist ) + c;
		err_vec = (data - est) .* weight;
		errs = sum(err_vec .^ 2);
	end
end
