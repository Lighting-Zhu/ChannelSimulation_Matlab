clear all;
EbN0dB = 0:10;%�趨�����
errLimit = 500;%��������ﵽ5ʱֹͣ����
N = 10^5;
M = 4;
k = log2(M);
divorder = 8; %�ּ�����
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

h = sqrt(0.5) * (randn(divorder,N) + 1j*randn(divorder,N) );

E = mean(abs(h).^2);



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
 n = 1/sqrt(2)*[randn(divorder,N) + 1j*randn(divorder,N)]; % white guassian noise, 0dB variance 
 y = ones(divorder,1)* s .* h + 10^(-EsN0dB(i)/20)*sqrt(divorder).*n; % additive white gaussian noise
 y = y./(h);
 
y = y.*abs(h).^2;
y = sum(y);
        %demode
        y_re = real(y).*sqrt(2);
        y_im = imag(y).*sqrt(2);
        
        
        demodx = zeros(1,N);
        demodx(y_re < 0) = 1;
        demody = zeros([1,N]);
        demody(y_im < 0) = 1;
        
        temp = sum(demodx ~= bitsx ) + sum(demody ~= bitsy);%���������
        errCount = errCount + temp;
        totalN = totalN + N;
    %end
    ber(i) = errCount/(totalN * k);
end
semilogy(EbN0dB, ber, 'o');%������ö�������

hold on;
grid on;

bertheory = berfading(EbN0dB,'psk',4,divorder);%���һ������Ϊ�ּ�����
semilogy(EbN0dB, bertheory, '-*');

axis([EbN0dB(1), EbN0dB(end), 10^-7 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              





        