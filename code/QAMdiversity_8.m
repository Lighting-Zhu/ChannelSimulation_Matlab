% Bit Error Rate for 16-QAM modulation using Gray modulation mapping

clear
N = 10^4; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol
divorder  = 8;%diversity order
% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-3 -1 1 3];
alphaIm = [-3 -1 1 3];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = 0:10; 
EsN0dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = [0:k-1];
%map = bitxor(ref,floor(ref/2));
%[tt, ind] = sort(map);                                
ind = [1 2 4 3];

for i = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    
    h = sqrt(0.5) * (randn(divorder,N) + 1j*randn(divorder,N) );
    E = mean(abs(h).^2);

    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';%ת��Ϊ4��
    % conversion from binary to decimalÿ��Ϊ2 1��������
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; 
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);%ȡ����1��2��
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));%��0123��Ϊ0132
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);%ӳ��Ϊ -3 -1 1 3
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + 1i*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*[randn(divorder,N) + 1j*randn(divorder,N)]; % white guassian noise, 0dB variance 
    y = ones(divorder,1)* s .* h + 10^(-EsN0dB(i)/20)*sqrt(divorder).*n; % additive white gaussian noise
    y = y./(h);
 
    y = y.*abs(h).^2;
    y = sum(y);
    y = y./sum(abs(h).^2) ;

    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet ���ն��źŹ�����[-3,3]
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(ipHatRe>max(alphaRe)) = max(alphaRe);
    ipHatRe(ipHatRe<min(alphaRe)) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(ipHatIm>max(alphaIm)) = max(alphaIm);
    ipHatIm(ipHatIm<min(alphaIm)) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based ӳ�䵽[1.5,4.5] ȡfloor
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string ʮ����ת��Ϊ������
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
    nBitErr(i) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);
%theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));
theoryBer = berfading(Eb_N0_dB,'qam',M,divorder);
semilogy(Eb_N0_dB,theoryBer,'-*');
hold on
semilogy(Eb_N0_dB,simBer,'o');
axis([Eb_N0_dB(1), Eb_N0_dB(end), 10^-7 1])

grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('16-QAM modulation')