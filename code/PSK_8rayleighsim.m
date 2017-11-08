clear all;
EbN0dB = 0:10;%设定信噪比
errLimit = 500;%误码个数达到5时停止仿真
N = 10^5;
sigma2 =  1./( 3 *  2 * 10.^(EbN0dB/10) );
ber = zeros(size(EbN0dB));
decibound = 1/8 : 1/4 : 7/8;
decibound = [decibound, -7/8:1/4:-1/8]; 
decibound = decibound * pi;

for i = 1:length(EbN0dB)
    errCount = 0;
    totalN = 0;
%     rand1 = wichhillrand(N);
%     rand2 = wichhillrand(N);
rand1 = rand(1,N);
rand2 = rand(1,N);
rand3 = rand(1,N);
h = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );

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
    complexsignal = exp(1j * phaseIQ );
    complexsignal = complexsignal .* h;
   
    I = real(complexsignal);
    Q = imag(complexsignal);
    
    
        %modulation
        %pass
        
        receiverI = I +sqrt(sigma2(i)) * randn(1,N) ;
        receiverQ = Q +sqrt(sigma2(i)) * randn(1,N) ;
        receiverIQ = (receiverI + 1j * receiverQ)./h;
        receiverI = real(receiverIQ);
        receiverQ = imag(receiverIQ);

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

semilogy(EbN0dB, ber , 'o');%纵轴采用对数坐标

hold on;
grid on;

bertheory = berfading(EbN0dB,'psk',8,1);
semilogy(EbN0dB, bertheory, '-*');

axis([EbN0dB(1), EbN0dB(end), 10^-2 1])
xlabel('Eb/N0 in dB')
ylabel('BER')
legend('simulation','theoretical')
                                                                                                                                                                                                                                              




