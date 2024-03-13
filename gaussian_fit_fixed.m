function y = gaussian_fit_fixed(b,x)

y = b(2) + b(1)*exp(-((x-0).^2)/2/30^2);