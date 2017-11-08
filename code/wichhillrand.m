function randarray = wichhillrand(a)
x0 = randi(200);
y0 = randi(200);
z0 = randi(200);
randarray = zeros(1,a);
x = zeros(1,a);
y = zeros(1,a);
z = zeros(1,a);
x(1) = mod(171 * x0, 30269);
y(1) = mod(170 * y0, 30307);
z(1) = mod(172 * z0, 30323);
randarray(1) = mod((x(1)/30269 + y(1)/30307 + z(1)/30323),1);
for i = 2:a
    x(i) = mod(171 * x(i-1),30269);
    y(i) = mod(170 * y(i-1),30307);
    z(i) = mod(172 * z(i-1),30323);
    randarray(i) = mod((x(i)/30269 + y(i)/30307 + z(i)/30323),1);
end
end
