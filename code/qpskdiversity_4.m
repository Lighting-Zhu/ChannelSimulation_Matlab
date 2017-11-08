clear all;
EbN0dB = 0:10;%设定信噪比
errLimit = 500;%误码个数达到5时停止仿真
N = 10^5;
M = 4;
k = log2(M);
divorder = 4; %分集阶数
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
E1 = mean(abs(h1).^2);
E2 = mean(abs(h2).^2);
E3 = mean(abs(h3).^2);
E4 = mean(abs(h4).^2);


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

 n1 = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; % white guassian noise, 0dB variance 
 y1 = s .* h1 + 10^(-EsN0dB(i)/20)*sqrt(divorder/E1)*n1; % additive white gaussian noise
 y1 = y1./(h1);
 
 n2 = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; % white guassian noise, 0dB variance 
 y2 = s .* h2 + 10^(-EsN0dB(i)/20)*sqrt(divorder/E2)*n2; % additive white gaussian noise
 y2 = y2./(h2);
 
  n3 = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; % white guassian noise, 0dB variance 
 y3 = s .* h3 + 10^(-EsN0dB(i)/20)*sqrt(divorder/E3)*n3; % additive white gaussian noise
 y3 = y3./(h3);
 
  n4 = 1/sqrt(2)*[randn(1,N) + 1j*randn(1,N)]; % white guassian noise, 0dB variance 
 y4 = s .* h4 + 10^(-EsN0dB(i)/20)*sqrt(divorder/E4)*n4; % additive white gaussian noise
 y4 = y4./(h4);
 %最大比合并，与理论曲线重合
y = y1.*abs(h1).^2/var(n1) + y2 .*abs(h2).^2./var(n2) ...
   + y3 .* abs(h3).^2./var(n3) + y4 .*abs(h4).^2./var(n4);
%等增益合并 BER更差
%y = y1.*abs(h1) + y2 .*abs(h2) + y3 .* abs(h3) + y4 .*abs(h4);
%选择性合并
% [~,index] = max(abs([h1; h2;h3;h4]));
% index = [0:4:4*N-1]+index;
% y = [y1;y2;y3;y4];
% y = y(:);
% y = y(index);
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
semilogy(EbN0dB, ber, '-o');%纵轴采用对数坐标

hold on;
grid on;

bertheory = berfading(EbN0dB,'psk',4,divorder);%最后一个参数为分集阶数
semilogy(EbN0dB, bertheory, '-*');

axis([EbN0dB(1), EbN0dB(end), 10^-7 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              



