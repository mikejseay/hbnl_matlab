function name = basename(input, suffix)
	name = '';
	remainder = input;
	while(any(remainder))
		[name, remainder] = strtok(remainder, '/');
	end
	if nargin == 2
		full_name = name;
		[name, remainder] = strtok(name, '.');
		while(any(remainder))
			[name, remainder] = strtok(remainder, '.');
		end
		if strcmp(name, suffix)
			name = full_name;
		else
			name = [];
		end
	end
