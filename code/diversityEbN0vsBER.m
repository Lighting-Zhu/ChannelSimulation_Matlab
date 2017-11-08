%%simulation using different diversity and modulation
%%Eb/N0 vs BER




%%
clear all;
EbN0dB = -10:10;%设定信噪比
N = 10^5;
k = log2(8);
sigma2 =  1./( 3 *  2 * 10.^(EbN0dB/10) );
ber = zeros(size(EbN0dB));
decibound = 1/8 : 1/4 : 7/8;
decibound = [decibound, -7/8:1/4:-1/8]; 
decibound = decibound * pi;
EsN0dB  = EbN0dB + 10*log10(k);

for i = 1:length(EsN0dB)
    errCount = 0;
    totalN = 0;
%     rand1 = wichhillrand(N);
%     rand2 = wichhillrand(N);
rand1 = rand(1,N);
rand2 = rand(1,N);
rand3 = rand(1,N);
h1 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h2 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h3 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h4 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );


    bitsx = zeros(1,N);
    bitsx(rand1>0.5) = 1;
    bitsy = zeros([1,N]);
    bitsy(rand2>0.5) = 1;
    bitsz = zeros([1,N]);
    bitsz(rand3>0.5) = 1;
    phaseIQ = zeros(1,N);
    phaseIQ(bitsx == 0 & bitsy == 0 & bitsz == 0) = 0;
    phaseIQ(bitsx == 0 & bitsy == 0 & bitsz == 1) = 1;
    phaseIQ(bitsx == 0 & bitsy == 1 & bitsz == 1) = 2;
    phaseIQ(bitsx == 0 & bitsy == 1 & bitsz == 0) = 3;
    phaseIQ(bitsx == 1 & bitsy == 1 & bitsz == 0) = 4;
    phaseIQ(bitsx == 1 & bitsy == 1 & bitsz == 1) = 5;
    phaseIQ(bitsx == 1 & bitsy == 0 & bitsz == 1) = 6;
    phaseIQ(bitsx == 1 & bitsy == 0 & bitsz == 0) = 7;
    phaseIQ = phaseIQ/4 * pi;
    s = exp(1j * phaseIQ );
     
    n1 = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance 
    % additive white gaussian noise 
    %分集每路信号Eb／No减小为分集数k的1/k即噪声乘以sqrt(k)
    y1 = s .* h1 + 10^(-EsN0dB(i)/20)*n1.*sqrt(2); 
    y1 = y1./(h1./abs(h1));
    
    n2 = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)];
    y2 = s .* h2 + 10^(-EsN0dB(i)/20)*n2.*sqrt(2); % additive white gaussian noise
    y2 = y2./(h2./abs(h2));
    
    y = (y1 .* abs(h1) + y2 .* abs(h2));
    
    
        %modulation
        %pass
       
  
        receiverI = real(y);
        receiverQ = imag(y);

        %demode
        dephase = angle(receiverI + 1j * receiverQ);
        demodxyz = zeros(1,N);
        demodxyz(dephase >= decibound(8) & dephase < decibound(1)) = 0;
        demodxyz(dephase >= decibound(1) & dephase < decibound(2)) = 1;
        demodxyz(dephase >= decibound(2) & dephase < decibound(3)) = 2;
        demodxyz(dephase >= decibound(3) & dephase < decibound(4)) = 3;
        demodxyz(dephase >= decibound(4) | dephase < decibound(5)) = 4;
        demodxyz(dephase >= decibound(5) & dephase < decibound(6)) = 5;
        demodxyz(dephase >= decibound(6) & dephase < decibound(7)) = 6;
        demodxyz(dephase >= decibound(7) & dephase < decibound(8)) = 7;

        xyz = zeros(3, length(demodxyz));
        xyz(:, demodxyz == 0) = [0; 0; 0] * ones(1,sum(demodxyz == 0));
        xyz(:, demodxyz == 1) = [0; 0; 1] * ones(1,sum(demodxyz == 1));
        xyz(:, demodxyz == 2) = [0; 1; 1] * ones(1,sum(demodxyz == 2));
        xyz(:, demodxyz == 3) = [0; 1; 0] * ones(1,sum(demodxyz == 3));
        xyz(:, demodxyz == 4) = [1; 1; 0] * ones(1,sum(demodxyz == 4));
        xyz(:, demodxyz == 5) = [1; 1; 1] * ones(1,sum(demodxyz == 5));
        xyz(:, demodxyz == 6) = [1; 0; 1] * ones(1,sum(demodxyz == 6));
        xyz(:, demodxyz == 7) = [1; 0; 0] * ones(1,sum(demodxyz == 7));

        
        
        demodx = xyz(1,:);
        demody = xyz(2,:);
        demodz = xyz(3,:);

        temp = sum(demodx ~= bitsx ) + sum(demody ~= bitsy) + sum((demodz ~= bitsz));
        errCount = errCount + temp;
        totalN = totalN + N * 3;
    %end
    ber(i) = errCount/totalN;
end

semilogy(EbN0dB, ber , '-o');%纵轴采用对数坐标

hold on;
grid on;

bertheory = berfading(EbN0dB,'psk',8,2);
semilogy(EbN0dB, bertheory, '-*');

axis([EbN0dB(1), EbN0dB(end), 10^-4 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              




