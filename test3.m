clear all;
clc;
close; 
%
%17206018 ChengLitao
% 1. The program should include the OFDM transmitter, 
% equivalent discrete-time channel, AWGN and OFDM demodulator. 
% 
% 2. For the equivalent discrete-time channel, 
% generate the channel filter coefficients h[n]Î¼n=0 
% as independent and identically distributed zero-mean complex Gaussian random variables, 
% with variance 1/2 for real and imaginary parts.
% 
% 3. For the OFDM demodulator you can assume that the channel coefficients h[n]Î¼n=0 are perfectly known.
% 
% 4. Test your MATLAB program with 64-QAM modulation on all data subcarriers and rc = 1/2 
% and obtain a plot of the bit error rate versus received Eb/N0.

SNR = 1:1:20; 
BER1 = zeros(1, length(SNR)); % bit error rate
subcarrierNum = 48; % the number of subcarrier with data
groupNum=2^12; 
bitsNum =6*subcarrierNum*groupNum; % data scale

% convolutional for rate 1/2 feedback
trellis = poly2trellis(5, [37 33], 37); 
tbdepth = 34; % Traceback depth

convertor = [32*ones(1, 2*bitsNum/6); 16*ones(1, 2*bitsNum/6); 8*ones(1, 2*bitsNum/6);4*ones(1, 2*bitsNum/6); 2*ones(1, 2*bitsNum/6); ones(1, 2*bitsNum/6)]; % bin2dec convertor

for snr = SNR
%1.CreatbitSignal
    bits_in = randi([0 1], bitsNum, 1);  
    
    %% OFDM Transmitter
%2convolutional code
    coded_bits = convenc(bits_in, trellis); 

%3.64-QAM modulations
    coded_bits = reshape(coded_bits, [6,  2*subcarrierNum*groupNum]); 
    coded_bits = sum(coded_bits.*convertor);
    coded_bits=reshape(coded_bits,[subcarrierNum,2*groupNum]);
    complex_signals = qammod(coded_bits, 64);
    
%4 ifft
    complex_signals = reshape(complex_signals, [subcarrierNum, 2*bitsNum/(subcarrierNum*6)]);
    % Perform 64 ifft operation
    % returns the inverse transform of each column of the matrix.
    %padding Y with trailing zeros to length n.
    complex_signals = ifft(complex_signals, 64);
    
%5 Add cyclic prefix
    cp_complex_signals = zeros(80, 2*bitsNum/(subcarrierNum*6)); 
    cp_complex_signals(17:end, :) = complex_signals; 
    cp_complex_signals(1:16, :) = complex_signals(49:64, :); 
    
    
    %% Equivalent discrete-time channel
%6.add AWGN
    % Generate the channel filter coefficients h[n]¦Ìn=0
    % as independent and identically distributed zero-mean complex Gaussian random variables, 
    % with variance 1/2 for real and imaginary parts.
    h = 1/(sqrt(0.5*randn+0.5*randn*1i)); 
    channel_rayleigh = h*cp_complex_signals; 
    noise_gaussian = awgn(channel_rayleigh, snr, 'measured'); 
    cp_complex_signals = h\noise_gaussian; 
    
    %% OFDM Receiver
%7remove CP
    complex_signals = cp_complex_signals(17:end, :); 
    
%8fft
    complex_signals = fft(complex_signals, 64);
    a=zeros(subcarrierNum,2*bitsNum/(subcarrierNum*6));
    a=complex_signals(1:48,:);
    
    
    
%9.64-QAM demodulation £¨optimum£©
    complex_signals = reshape(a, [2*bitsNum/6, 1]); 
    coded_bits = qamdemod(a, 64); 
    coded_bits = dec2bin(coded_bits, 6)=='1'; 
    
%10remove convolutional code
    coded_bits = reshape(coded_bits', [bitsNum*2 ,1]); 
    bit_out = vitdec(coded_bits, trellis, tbdepth, 'trunc', 'hard'); 
    
%11calculate BER
    [number, ratio] = biterr(bit_out, bits_in); 
    BER1(snr) = ratio; 
end


%%without zeros channels simulation

SNR = 1:1:20; 
BER2 = zeros(1, length(SNR)); % bit error rate
subcarrierNum = 64; 
groupNum=2^12;
bitsNum =6*subcarrierNum*groupNum; % data scale

% convolutional for rate 1/2 feedback
trellis = poly2trellis(5, [37 33], 37); 
tbdepth = 34; % Traceback depth

convertor = [32*ones(1, 2*bitsNum/6); 16*ones(1, 2*bitsNum/6); 8*ones(1, 2*bitsNum/6);4*ones(1, 2*bitsNum/6); 2*ones(1, 2*bitsNum/6); ones(1, 2*bitsNum/6)]; % bin2dec convertor

for snr = SNR
%1.CreatbitSignal
    bits_in = randi([0 1], bitsNum, 1); 
    
    %% OFDM Transmitter
%2convolutional code
    %for i=1:1:subcarrierNum*groupNum
    coded_bits = convenc(bits_in, trellis); 
    %end

%3.64-QAM modulations
    coded_bits = reshape(coded_bits, [6, 2*bitsNum/6]); 
    coded_bits = sum(coded_bits.*convertor);
    complex_signals = qammod(coded_bits, 64); 
    
%4 ifft
    complex_signals = reshape(complex_signals, [subcarrierNum, 2*bitsNum/(subcarrierNum*6)]); % reshape to a matrix who has 2^i row
    
    % Perform 2^i-point ifft operation
    complex_signals = ifft(complex_signals, subcarrierNum); 
    
%5 Add cyclic prefix
    cp_complex_signals = zeros(80, 2*bitsNum/(subcarrierNum*6)); 
    cp_complex_signals(17:end, :) = complex_signals; 
    cp_complex_signals(1:16, :) = complex_signals(49:64, :); 
    
    
    %% Equivalent discrete-time channel
%6.add AWGN
    % Generate the channel filter coefficients h[n]¦Ìn=0
    % as independent and identically distributed zero-mean complex Gaussian random variables, 
    % with variance 1/2 for real and imaginary parts.
    h = 1/(sqrt(0.5*randn+0.5*randn*1i)); 
    channel_rayleigh = h*cp_complex_signals; 
    noise_gaussian = awgn(channel_rayleigh, snr, 'measured'); 
    cp_complex_signals = h\noise_gaussian; 
    
    %% OFDM Receiver
%7remove CP
    complex_signals = cp_complex_signals(17:end, :); 
    
%8fft
    complex_signals = fft(complex_signals, subcarrierNum); 
    
%9.64-QAM demodulation £¨optimum£©
    complex_signals = reshape(complex_signals, [2*bitsNum/6, 1]); 
    coded_bits = qamdemod(complex_signals, 64); 
    coded_bits = dec2bin(coded_bits, 6)=='1'; 
    
%10remove convolutional code
    coded_bits = reshape(coded_bits', [bitsNum*2 ,1]); 
    dec_bits = vitdec(coded_bits, trellis, tbdepth, 'trunc', 'hard'); 
    
%11calculate BER
    [number, ratio] = biterr(dec_bits, bits_in); 
    BER2(snr) = ratio; 
end


% plot
figure(1);
semilogy(SNR+10.*log10(6), BER2,'-b.','MarkerSize',10); 
hold on;
semilogy(SNR+10.*log10(6), BER1,'-r.','MarkerSize',10); 
title("OFDM system BER performances")
xlabel("Eb/N0 (dB)")
ylabel("BitErrorRate")
legend("without zeros channels","64-QAM BER")
grid on;
