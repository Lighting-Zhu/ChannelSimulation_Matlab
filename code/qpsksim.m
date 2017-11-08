clear all;
EbN0dB = -10:10;%设定信噪比
errLimit = 500;%误码个数达到5时停止仿真
N = 1000000;
sigma2 =  1./( 2 *10.^(EbN0dB/10) );
ber = zeros(size(EbN0dB));


for i = 1:length(EbN0dB)
    errCount = 0;
    totalN = 0;
%     rand1 = wichhillrand(N);
%     rand2 = wichhillrand(N);
rand1 = rand(1,N);
rand2 = rand(1,N);
    bitsx = zeros([1,N]);
    bitsx(rand1>0.5) = 1;
    bitsy = zeros([1,N]);
    bitsy(rand2>0.5) = 1;
    symsx = 1 - 2 * bitsx;
    symsy = 1 - 2 * bitsy;
    
        %modulation
        %pass
        
        receiverx = symsx +sqrt(sigma2(i)) * randn(1,N);
        receivery = symsy +sqrt(sigma2(i)) * randn(1,N);

        %demode
        demodx = zeros(1,N);
        demodx(receiverx < 0) = 1;
        demody = zeros([1,N]);
        demody(receivery < 0) = 1;
        
        temp = sum(demodx ~= bitsx ) + sum(demody ~= bitsy);%错误比特数
        errCount = errCount + temp;
        totalN = totalN + N;
    %end
    ber(i) = errCount/totalN;
end
semilogy(EbN0dB, ber /2, '-o');%纵轴采用对数坐标

hold on;
grid on;

bertheory = berawgn(EbN0dB,'psk',4,'nodiff');
semilogy(EbN0dB, bertheory, '-+');

axis([EbN0dB(1), EbN0dB(end), 10^-7 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              





        