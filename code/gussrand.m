function x = gussrand(sigma, num)
%��Ҫ���������������ʹ�ò�ͬ������ sigma: std

r = rayleigh(sigma, num);
u2 = wichhillrand(num);
x = r .* cos(2 * pi * u2);


end