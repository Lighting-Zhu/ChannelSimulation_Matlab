clear all;
EbN0dB = -10:10;%设定信噪比
errLimit = 500;%误码个数达到5时停止仿真
N = 100000;
sigma2 = 1./( 2 * 10.^(EbN0dB/10));
ber = zeros(size(EbN0dB));
rand1 = wichhillrand(N);
gussrand1 = gussrand(1,2*N);

for i = 1:length(EbN0dB)
    errCount = 0;
    totalN = 0;
    %while errCount<errLimit
    bits = zeros([1,N]);
        bits(rand1>0.5) = 1;
        
        syms = 1 - 2 * bits;

        %modulation
        %pass
        
        receiver = syms +sqrt(sigma2(i)) * gussrand1(N + 1:2 * N);
        
        %demode
        demod = zeros(1,N);
        demod(receiver < 0) = 1;
        
        temp = sum(abs(bits-demod));
        errCount = errCount + temp;
        totalN = totalN + N;
    %end
    ber(i) = errCount/totalN;
end
semilogy(EbN0dB, ber, '-o');%纵轴采用对数坐标

hold on;
grid on;

bertheory = berawgn(EbN0dB,'psk',2,'nodiff');
semilogy(EbN0dB, bertheory, '-*');

axis([EbN0dB(1), EbN0dB(end), 10^-10 1])
xlabel('SNR in dB')
ylabel('BER')
legend('simulation','theoretical')






        