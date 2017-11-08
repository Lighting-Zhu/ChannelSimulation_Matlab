function x = gussrand(sigma, num)
%需要产生两组随机数，使用不同的种子 sigma: std

r = rayleigh(sigma, num);
u2 = wichhillrand(num);
x = r .* cos(2 * pi * u2);


end