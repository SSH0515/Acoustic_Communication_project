clc
clear all
close all

bits=[1 0 1 0 1 1 1 0 1 0 1 0 1 1 1 0]; % 예시 데이터


%% CRC generator
poly = [1 0 1 1];
reshaped_btis = reshape(bits,8,length(bits)/8);
reshaped_transposed_bits = reshaped_btis.';
add_zero_bits = [reshaped_transposed_bits zeros(height(reshaped_transposed_bits),3)];

for i = 1:height(add_zero_bits)
    remainder_Tx = [reshaped_transposed_bits(i,:) 0 0 0];
    for j = 1:length(reshaped_transposed_bits)
        if remainder_Tx(j) == 1
            remainder_Tx(j:j+length(poly)-1) = xor(remainder_Tx(j:j+length(poly)-1),poly);
        end
    end
    CRC = remainder_Tx(end-2:end);
    add_zero_bits(i,:) = [reshaped_transposed_bits(i,:) CRC];
end

b = add_zero_bits.';
bits_w_CRC = b(:);

%% CRC generator
% poly = [1 0 1 1]; % X^3 + X^2 + 1
% 
% remainder_Tx = [bits, 0 0 0];
% for i = 1:length(bits)
%     if remainder_Tx(i)==1
%         % XOR 연산 수행 (다항식 나눗셈)
%         remainder_Tx(i:i+length(poly)-1) = xor(remainder_Tx(i:i+length(poly)-1),poly);
%     end
% end
% 
% CRC = remainder_Tx(end-2:end); % 나머지 부분
% 
% bits_w_CRC = [bits CRC];

%% Transmitter
k = 11; % Number of valid bits in 1 coded bits
c = 4; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
P = [1 0 1 1;1 1 1 0;1 1 0 1;1 1 0 0;0 1 1 1;1 0 1 1;1 1 1 1;0 1 0 1;1 0 1 0;0 1 1 0;1 0 0 1]; % Parity generator matrix
G = [eye(k),P]; % code generator matrix

reshaped_bits = reshape(bits_w_CRC,[k,length(bits_w_CRC)/k]);
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

%% CRC check
crc_bits = reshape(HD_decoded_bits,11,length(HD_decoded_bits)/11).';
CRC_check_res = [];

for i = 1:height(crc_bits)
    remainder_Rx = crc_bits(i,:);
    for j = 1:size(crc_bits,2)-length(poly)+1
        if remainder_Rx(j) == 1
            remainder_Rx(j:j+length(poly)-1) = xor(remainder_Rx(j:j+length(poly)-1),poly);
        end
    end
    CRC_check_res = [CRC_check_res; remainder_Rx(end-2:end)];
end

%% CRC Check
%HD_decoded_bits(4)=1; % 에러난 상황 가정
% remainder_Rx=HD_decoded_bits;
% for i = 1:length(HD_decoded_bits)-length(poly)+1
%     if remainder_Rx(i)==1
%         % XOR 연산 수행 (다항식 나눗셈)
%         remainder_Rx(i:i+length(poly)-1) = xor(remainder_Rx(i:i+length(poly)-1),poly);
%     end
% end
% 
% % 계산된 나머지
% CRC_check = remainder_Rx(end-2:end); % 나머지
% 
% % 오류 검증
% if all(CRC_check ==0)
%     disp('NO ERROR');
% else
%     disp('ERROR');
% end