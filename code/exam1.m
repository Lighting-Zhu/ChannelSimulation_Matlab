%%
rand1 = wichhillrand(10^4);
rand2 = wichhillrand(10^4);
figure(1)
hist(rand1,20);
figure(2)
hist(rand2,20);

%%
rand3 = rayleigh(sqrt(0.5), 10^4);
rand4 = rayleigh(sqrt(1), 10^4);
rand5 = rayleigh(sqrt(2), 10^4);
figure(3)
hist(rand3,20);
figure(4)
hist(rand4,20);
figure(5)
hist(rand5,20);

%%
rand6 = gussrand(1,10^5);
rand7 = gussrand(3,10^5);
figure(6)
hist(rand6,20);
figure(7)
hist(rand7,20);








