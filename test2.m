clear all;
clc;
close; 
%%
%17206018 ChengLitao
%tasks:
%
%
%

SNR = 1:1:20; 
BER = zeros(1, length(SNR)); 
trellis = poly2trellis(5, [37 33], 37); 
tbdepth = 34; % Traceback depth for Viterbi decoder
groupNum=6*100;
subcarrierNum=64;
convertor = [32*ones(1, groupNum*subcarrierNum); 16*ones(1, groupNum*subcarrierNum); 8*ones(1, groupNum*subcarrierNum); ...
4*ones(1, groupNum*subcarrierNum); 2*ones(1, groupNum*subcarrierNum); ones(1, groupNum*subcarrierNum)]; 

for snr = 1:1:20
%1.CreatbitSignal
bits_in=rand(groupNum*subcarrierNum*6/2,1)<0.5;

%2convolutional code
enc_bits = convenc(bits_in, trellis);

%3.64-QAM modulations
enc_bits=reshape(enc_bits,[6,groupNum*subcarrierNum]);
enc_bits = sum(enc_bits.*convertor);
complex_signals = qammod(enc_bits, 64); 

%4 ifft
complex_signals=reshape(complex_signals,[subcarrierNum,groupNum]);
complex_signals=ifft(complex_signals,subcarrierNum);

%5 add CP 
cp_complex_signals = zeros(80,groupNum); 
cp_complex_signals(17:end, :) = complex_signals; 
cp_complex_signals(1:16, :) = complex_signals(49:64, :);
 
%6.add AWGN
% Generate the channel filter coefficients h[n]w n=0
% as independent and identically distributed zero-mean complex Gaussian random variables,
% with variance 1/2 for real and imaginary parts.
h = 1/(sqrt(0.5*randn+0.5*randn*1i));
channel_rayleigh = h*cp_complex_signals;
noise_gaussian = awgn(channel_rayleigh, snr, 'measured');
cp_complex_signals = h\noise_gaussian;  %inv(h)*noise_gaussian


%7remove CP
complex_signals = cp_complex_signals(17:end, :); 

%8fft
complex_signals = fft(complex_signals, subcarrierNum);

%9.64-QAM demodulation £¨optimum£©
enc_bits = qamdemod(complex_signals, 64); 
enc_bits = dec2bin(enc_bits, 6)=='1'; 

%10remove convolutional code
enc_bits = reshape(enc_bits', [groupNum*subcarrierNum*6,1]); 
bit_out = vitdec(enc_bits, trellis, tbdepth, 'trunc', 'hard');

%calculate BER
[number, ratio] = biterr(bit_out, bits_in); 
BER(snr) = ratio;
end

% plot
figure(1);
semilogy(SNR/6, BER,'r'); 
title("OFDM system bit error rate performances")
xlabel("SNR (Eb/N0)")
ylabel("BitErrorRate")
grid on;

