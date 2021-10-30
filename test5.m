clear all;
clc;
close; 
%
%17206018 ChengLitao


%tasks:
% arguments and constants

SNR = 1:1:20; 
BER2 = zeros(1, length(SNR)); % bit error rate
subcarrierNum = 64; 
groupNum=1000;
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

BER1=[0.368136088053385	0.352842542860243	0.335761176215278	0.318423800998264	0.297899034288194	0.274937947591146	0.248069763183594	0.214982774522569	0.174608866373698	0.128533257378472	0.0814666748046875	0.0434705946180556	0.0195812649197049	0.00724199083116319	0.00258721245659722	0.000714619954427083	0.000189887152777778	4.15378146701389e-05	3.39084201388889e-06	0];

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


