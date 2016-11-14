n_samps = 384;
rate = 256;
freq_lims = [2 32];
padratio = 2;
[scales_vec, cycles] = findCWTscales(n_samps, rate, freq_lims, padratio);
cycles_vec = logspace(log10(cycles(2)),log10(cycles(1)),length(scales_vec));
n_scales = length(scales_vec);

%cycle_type = 'fixed';
cycle_type = 'variable';

figure;
for scale_ind=1:n_scales

scale = scales_vec(scale_ind);
n_cycles = cycles_vec(scale_ind);

switch cycle_type
    
    case 'fixed'
        width = 4;
        N = width * scale; %length in samples
        wo = 2 * pi;
        t = linspace(-width, width, N); %analytical time
        cycle_const = (n_cycles+1) ./ (2 * pi);
        gauss = exp(-(t .^ 2)./2);
        ww = sqrt(1/(scale)) .* exp(1i * wo .* t) .* gauss;
        wwf = @(x) sqrt(1/(scale)) .* exp(1i * wo .* x) .* exp(-(x .^ 2)./2);
    
    case 'variable'
        width = round( 2 .* n_cycles );
        N = width * scale; %length in samples
        wo = 2 * pi;
        t = linspace(-width, width, N); %analytical time
        cycle_const = (n_cycles+1) ./ (2 * pi);
        gauss = exp(-(t .^ 2)./(2.*cycle_const.^2));
        ww = sqrt(1/(scale.*cycle_const)) .* exp(1i * wo .* t) .* gauss;
        wwf = @(x) sqrt(1/(scale.*cycle_const)) .* exp(1i * wo .* x) .* ...
            exp(-(x .^ 2)./(2.*cycle_const.^2));
end

integ(scale_ind) = integral(wwf, -Inf, Inf);
integ2(scale_ind) = sum(ww);

subplot(4,4,scale_ind);
plot(real(ww)); title(num2str(n_cycles));
axis([1 N -0.25 0.25]);

end

figure;
subplot(211); plot(real(integ));
subplot(212); plot(real(integ2));