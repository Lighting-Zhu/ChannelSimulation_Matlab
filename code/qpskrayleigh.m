clear all;
EbN0dB = 0:10;%设定信噪比
errLimit = 500;%误码个数达到5时停止仿真
N = 10^5;
M = 4;
k = log2(M);
divorder = 1; %分集阶数
sigma2 =  1./( 2 *10.^(EbN0dB/10) );
ber = zeros(size(EbN0dB));
EsN0dB  = EbN0dB + 10*log10(k);
%EsN0dB = EsN0dB - 10*log10(divorder);


for i = 1:length(EbN0dB)
    errCount = 0;
    totalN = 0;
%     rand1 = wichhillrand(N);
%     rand2 = wichhillrand(N);
rand1 = rand(1,N);
rand2 = rand(1,N);
h1 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h2 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h3 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h4 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );

    bitsx = zeros([1,N]);
    bitsx(rand1>0.5) = 1;
    bitsy = zeros([1,N]);
    bitsy(rand2>0.5) = 1;
    symsx = 1 - 2 * bitsx;
    symsy = 1 - 2 * bitsy;
    s = symsx + 1j * symsy;
    s = s./sqrt(2); % normalization
    
        %modulation
        %pass

 n1 = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)].*sqrt(divorder); % white guassian noise, 0dB variance 
 y1 = s .* h1 + 10^(-EsN0dB(i)/20)*n1; % additive white gaussian noise
 y1 = y1./(h1)*sqrt(2);
 
 
y = y1;


        %demode
        y_re = real(y).*sqrt(2);
        y_im = imag(y).*sqrt(2);
        
        
        demodx = zeros(1,N);
        demodx(y_re < 0) = 1;
        demody = zeros([1,N]);
        demody(y_im < 0) = 1;
        
        temp = sum(demodx ~= bitsx ) + sum(demody ~= bitsy);%错误比特数
        errCount = errCount + temp;
        totalN = totalN + N;
    %end
    ber(i) = errCount/(totalN * k);
end

bertheory = berfading(EbN0dB,'psk',4,divorder);%最后一个参数为分集阶数
semilogy(EbN0dB, bertheory, '-*');
hold on;
grid on;

semilogy(EbN0dB, ber, 'o');%纵轴采用对数坐标

axis([EbN0dB(1), EbN0dB(end), 10^-7 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              





        