clc
clear all
close all

A = zeros(1,11);
B = zeros(1,11);
A(1,6)=1;
B(1,11)=1;

bits = [A, B];

%% Trnasmitter
k = 11; % Number of valid bits in 1 coded bits
c = 4; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
P = [1 0 1 1;1 1 1 0;1 1 0 1;1 1 0 0;0 1 1 1;1 0 1 1;1 1 1 1;0 1 0 1;1 0 1 0;0 1 1 0;1 0 0 1]; % Parity generator matrix
G = [eye(k),P]; % code generator matrix

reshaped_bits = reshape(bits,[k,length(bits)/k]);
haming_coded_bits = mod(transpose(reshaped_bits)*G,2);

transpose_haming_coded_bits = haming_coded_bits.';
channel_coded_bits = transpose_haming_coded_bits(:);

transmit_symbols = 2 *channel_coded_bits - 1; % BPSK Modulation

%% Receiver (hard decesion w/o noise)
demodulated_bits = (sign(transmit_symbols)+1)/2; % BPSK Demodulation

H =[P' eye(c)]; % Decoding matrix
reshaped_demodulated_bits = reshape(demodulated_bits, [n,length(demodulated_bits)/n]).';
syndrome_matrix = mod(reshaped_demodulated_bits*H.',2);

HD_decoded_bits = [];
for ii = 1:height(syndrome_matrix)
    find = 0; % index를 기록하기 위한 변수
    H_transpose = H.';
    for jj = 1:n
        if syndrome_matrix(ii,:) == H_transpose(jj,:) % syndrome을 이용하여 몇 번째 비트에서 error가 발생했는지 확인
            find = jj;
        end
    end
    if find == 0
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:11)]; %error가 발생하지 않은 경우, 받은 비트 그대로 사용
    else
        reshaped_demodulated_bits(ii,find) = mod((reshaped_demodulated_bits(ii,find)+1),2); % 에러가 발생한 곳 correction
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:11)];
    end
end

HD_decoded_bits