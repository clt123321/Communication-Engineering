%% Assignment Requirements
% 1. The program should include the OFDM transmitter, 
% equivalent discrete-time channel, AWGN and OFDM demodulator.
% 2. For the equivalent discrete-time channel, 
% generate the channel filter coefficients h[n]��n=0 as independent and identically distributed zero-mean complex Gaussian random variables, 
% with variance 1/2 for real and imaginary parts.
% 3. For the OFDM demodulator you can assume that the channel coefficients h[n]��n=0 are perfectly known.
% 4. Test your MATLAB program with 64-QAM modulation on all data subcarriers and rc = 1/2 
% and obtain a plot of the bit error rate versus received Eb/N0.


%% Simulation
clear; clc; 
t1 = clock; 

% arguments
Nsp = 52; %ϵͳ���ز�������������ֱ���ز���
Nfft = 64; %FFT����
Ncp = 16; %ѭ��ǰ׺����
Ns = Nfft+Ncp; %һ������OFDM���ų���
noc = 53; %����ֱ���ز����ܵ����ز���
Nd = 6; %ÿ֡������OFDM��������������ѵ�����ţ�
M1 = 4; %QPSK����
M2 = 16; %16-QAM����
sr = 250000; %OFDM��������
EbN0 = 0:2:30; %��һ�������
Nfrm = 10000; %ÿ��������µķ���֡��
ts = 1/sr/Ns; %OFDM���ų���ʱ����
t = 0:ts:(Ns*(Nd+1)*Nfrm-1)*ts; %����ʱ��
fd = 100; %��������Ƶ��
h = rayleigh(fd, t); %���ɵ���Rayleigh˥���ŵ�

%ѵ������Ƶ�����ݣ�����802.11a�еĳ�ѵ����������
Preamble = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ... 
            1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1]; 
Preamble1 = zeros(1,Nfft); 
Preamble1(2:27) = Preamble(27:end); %ѵ���������ź������??? 
Preamble1(39:end) = Preamble(1:26); 
preamble1 = ifft(Preamble1); %ѵ������ʱ������??? 
preamble1 = [preamble1(Nfft-Ncp+1:end) preamble1]; %����ѭ��ǰ׺

% preallocation for speed
ber1 = zeros(1, length(EbN0)); 
ber2 = zeros(1, length(EbN0)); 
ber3 = zeros(1, length(EbN0)); 
ber4 = zeros(1, length(EbN0)); 

for ii = 1:length(EbN0)
    %%%���������%%%
    msg1 = randsrc(Nsp, Nd*Nfrm, [0:M1-1]); %QPSK��Ϣ����?? 
    msg2 = randsrc(Nsp, Nd*Nfrm, [0:M2-1]); %16-QAM��Ϣ����?? 
    data1 = pskmod(msg1, M1, pi/4); %QPSK����?? 
    data2 = qammod(msg2, M2)/sqrt(10); %16-QAM���Ʋ���һ��
    
    data3 = zeros(Nfft, Nd*Nfrm); %����FFTҪ�󣬶���������?????????? 
    data4 = zeros(Nfft, Nd*Nfrm); 
    
    data3(2:27, :) = data1(27:end, :); %���ź������
    data3(39:end, :) = data1(1:26, :); 
    
    data4(2:27, :) = data2(27:end, :); %���ź������
    data4(39:end, :) = data2(1:26, :); 
    
    clear data1 data2; %�������Ҫ����ʱ����? 
    
    % 3. IFFT
    data3 = ifft(data3); 
    data4 = ifft(data4); 
    
    % 4. Add cyclic prefix
    data3 = [data3(Nfft-Ncp+1:end, :); data3]; 
    data4 = [data4(Nfft-Ncp+1:end, :); data4]; 
    
    spow1 = norm(data3, 'fro').^2/(Nsp*Nd*Nfrm); %�������ݷ�������??????????? 
    spow2 = norm(data4, 'fro').^2/(Nsp*Nd*Nfrm); % norm(X,'fro') returns the Frobenius norm of X.
    
    data5 = zeros(Ns, (Nd+1)*Nfrm); %����ѵ������??? 
    data6 = data5; 
    
    for indx = 1:Nfrm 
        data5(:, (indx-1)*(Nd+1)+1) = preamble1.'; 
        data5(:, (indx-1)*(Nd+1)+2:indx*(Nd+1)) = data3(:, (indx-1)*Nd+1:indx*Nd); 
        
        data6(:, (indx-1)*(Nd+1)+1) = preamble1.'; 
        data6(:, (indx-1)*(Nd+1)+2:indx*(Nd+1)) = data4(:, (indx-1)*Nd+1:indx*Nd); 
    end 
    
    clear data3 data4; 
    
    data5 = reshape(data5, 1, Ns*(Nd+1)*Nfrm); %�����任  
    data6 = reshape(data6, 1, Ns*(Nd+1)*Nfrm); 
    
    sigma1 = sqrt(1/2*spow1/log2(M1)*10.^(-EbN0(ii)/10)); %����EbN0����������׼��????
    sigma2 = sqrt(1/2*spow2/log2(M2)*10.^(-EbN0(ii)/10)); 
    
    % preallocation for speed
    neb1 = zeros(1, Nfrm); 
    neb2 = zeros(1, Nfrm); 
    neb3 = zeros(1, Nfrm); 
    neb4 = zeros(1, Nfrm); 
    
    for indx = 1:Nfrm 
        dd1 = data5((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1)); %��ǰ֡�ķ�������????? 
        dd2 = data6((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1)); 
        
        hh = h((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1)); %��ǰ֡��Ӧ���ŵ�����??  
        
        %�ź�ͨ������Rayleigh˥���ŵ����������˹������???????????????? 
        r1 = hh.*dd1+sigma1*(randn(1, length(dd1))+j*randn(1, length(dd1))); 
        r2 = hh.*dd2+sigma2*(randn(1, length(dd2))+j*randn(1, length(dd2))); 
        
        r1 = reshape(r1, Ns, Nd+1); %�����任
        r2 = reshape(r2, Ns, Nd+1); 
        
        r1 = r1(Ncp+1:end, :); %�Ƴ�ѭ��ǰ׺
        r2 = r2(Ncp+1:end, :);     
        
        
        %%%�����ŵ�����%%%
        hh = reshape(hh, Ns, Nd+1); %�ŵ�������������????????? 
        hh = hh(Ncp+1:end, :); 
        x1 = r1(:, 2:end)./hh(:, 2:end); %�����ŵ����ƽ��
        x2 = r2(:, 2:end)./hh(:, 2:end); 
        
        % 3. FFT
        x1 = fft(x1); 
        x2 = fft(x2); 
        
        x1 = [x1(39:end, :);x1(2:27, :)]; %��������??????? 
        x2 = [x2(39:end, :);x2(2:27, :)]; 
        
        % 4. Demodulation
        x1 = pskdemod(x1, M1, pi/4); 
        x2 = qamdemod(x2*sqrt(10), M2); 
        [neb1(indx), ~] = biterr(x1, msg1(:, (indx-1)*Nd+1:indx*Nd), log2(M1)); %ͳ��һ֡�еĴ����?????? 
        [neb2(indx), ~] = biterr(x2, msg2(:, (indx-1)*Nd+1:indx*Nd), log2(M2)); 
        
        %%%����ѵ�����Ž��е��ŵ�����%%%
        
        R1 = fft(r1); %fft����? 
        R2 = fft(r2); 
        
        R1 = [R1(39:end, :);R1(2:27, :)]; %��������?????? 
        R2 = [R2(39:end, :);R2(2:27, :)]; 
        
        HH1 = (Preamble.')./R1(:, 1); %�ŵ����� ??? 
        HH2 = (Preamble.')./R2(:, 1);                    
        
        HH1 = HH1*ones(1, Nd); %�����ŵ����ƽ�������ŵ�����?? 
        HH2 = HH2*ones(1, Nd); 
        
        x3 = R1(:, 2:end).*HH1;                         
        x4 = R2(:, 2:end).*HH2; 
        
        x3 = pskdemod(x3, M1, pi/4); %���ݽ��???? 
        x4 = qamdemod(x4.*sqrt(10), M2);         
        
        [neb3(indx), ~] = biterr(x3, msg1(:, (indx-1)*Nd+1:indx*Nd), log2(M1)); %ͳ��һ֡�еĴ����?????? 
        [neb4(indx), ~] = biterr(x4, msg2(:, (indx-1)*Nd+1:indx*Nd), log2(M2)); 
    end
    
    ber1(ii) = sum(neb1)/(Nsp*log2(M1)*Nd*Nfrm); %�����ŵ����Ƶ��������????? 
    ber2(ii) = sum(neb2)/(Nsp*log2(M2)*Nd*Nfrm); 
    
    ber3(ii) = sum(neb3)/(Nsp*log2(M1)*Nd*Nfrm); %����ѵ�������ŵ����Ƶ��������??????? 
    ber4(ii) = sum(neb4)/(Nsp*log2(M2)*Nd*Nfrm);          
end 

semilogy(EbN0, ber1, '-ro', EbN0, ber3, '-rv', EbN0, ber2, '-g*', EbN0, ber4, '-gp')
grid on
title("OFDMϵͳ�����������")
legend("QPSK�����ŵ�����", "QPSKѵ�������ŵ�����", "16-QAM�����ŵ�����", "16-QAMѵ�������ŵ�����")
xlabel("����ȣ�EbN0��")
ylabel("�������")

t2 = clock; 
etime(t2, t1)
