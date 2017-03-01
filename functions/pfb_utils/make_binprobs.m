function props = make_binprobs(n_bins)

% start between 0.1 and 0.9
% each change will be either positive or negative
% the change in probability will be sampled uniformly from [.2 .4]
% number of upward and downward changes in probability will be kept equal
% number of probabilities > and < 0.5 will be kept roughly equal

props = zeros(n_bins, 1);
props(1) = 0.1 + 0.8*rand;

while sum(sign(diff(props))) ~= 0 || abs(sum(sign(props - 0.5))) > 1

props(1) = 0.1 + 0.8*rand; % start somewhere in [0.1 0.9]

for bin=2:n_bins
    
    while (props(bin) < 0.1) || (props(bin) > 0.9) || (props(bin) > 0.45 && props(bin) < 0.55)
        direction = (randperm(2,1) - 1.5) * 2; % returns 1 or -1
        magnitude = 0.2 + rand*0.2; %returns [0.2 0.4]

        props(bin) = props(bin-1) + direction*magnitude;
    end
    
end

end

end