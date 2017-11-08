% Bit Error Rate for 16-QAM modulation using Gray modulation mapping

clear
N = 10^4; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol

% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-3 -1 1 3];
alphaIm = [-3 -1 1 3];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = [0:10]; 
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = [0:k-1];
%map = bitxor(ref,floor(ref/2));
%[tt, ind] = sort(map);                                
ind = [1 2 4 3];

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';%转化为4列
    % conversion from binary to decimal每行为2 1的列向量
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; 
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);%取出第1、2列
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));%将0123变为0132
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);%映射为 -3 -1 1 3
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + 1i*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance 
    
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    % ------------
    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet 将收端信号规整到[-3,3]
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(ipHatRe>max(alphaRe)) = max(alphaRe);
    ipHatRe(ipHatRe<min(alphaRe)) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(ipHatIm>max(alphaIm)) = max(alphaIm);
    ipHatIm(ipHatIm<min(alphaIm)) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based 映射到[1.5,4.5] 取floor
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string 十进制转化为二进制
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % counting errors for real and imaginary
    nBitErr(ii) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);
%theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));
theoryBer = berawgn(Eb_N0_dB,'qam',M);

semilogy(Eb_N0_dB,theoryBer,'-*');
hold on
semilogy(Eb_N0_dB,simBer,'o');

axis([Eb_N0_dB(1), Eb_N0_dB(end), 10^-7 1])

grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-QAM modulation')