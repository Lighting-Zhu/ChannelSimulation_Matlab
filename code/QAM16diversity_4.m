% Bit Error Rate for 16-QAM modulation using Gray modulation mapping

clear
N = 10^4; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol
divorder  = 4;
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
    
h1 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h2 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h3 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
h4 = sqrt(0.5) * (randn(1,N) + 1j*randn(1,N) );
E1 = mean(abs(h1).^2);
E2 = mean(abs(h2).^2);
E3 = mean(abs(h3).^2);
E4 = mean(abs(h4).^2);

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
%  y = y1.*abs(h1).^2/var(n1) + y2 .*abs(h2).^2./var(n2) + y3 .* abs(h3).^2./var(n3) + y4 .*abs(h4).^2./var(n4);
%  y = y./(abs(h1).^2/var(n1) + abs(h2).^2/var(n2) + abs(h3).^2/var(n3)  +abs(h4).^2./var(n4));
%等增益合并
 y = y1.*abs(h1) + y2 .*abs(h2) + y3 .* abs(h3) + y4 .*abs(h4);
 y = y./(abs(h1)+abs(h2)+abs(h3)+abs(h4));
 %选择性合并
[~,index] = max(abs([h1; h2;h3;h4]));
index = [0:4:4*N-1]+index;
y = [y1;y2;y3;y4];
y = y(:);
y = y(index);


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
    nBitErr(i) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);

% theoryBer = berfading(Eb_N0_dB,'qam',M,divorder);
% semilogy(Eb_N0_dB,theoryBer,'-*');
hold on
semilogy(Eb_N0_dB,simBer,'-o');
axis([Eb_N0_dB(1), Eb_N0_dB(end), 10^-7 1])

grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('16-QAM modulation')