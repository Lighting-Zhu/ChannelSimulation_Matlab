function r = rayleigh(sigma, num)
%sigma: std
a = wichhillrand(num);
r = sqrt(-2 * sigma^2 * log(a));

end