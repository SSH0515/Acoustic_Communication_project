clc
clear all
close all

%% Preamble
omega = 10;
mu =0.1;
Tp = 100;
tp = (1:Tp).';
preamble = cos(omega*tp+mu*tp.^2/2);

%% Source recording
img_Lena = imread('Lena_color.png'); % Load image
bits_Lena = imagetoBits(img_Lena);

%% Add CRC
bits_Lena_with_CRC = ADD_CRC(bits_Lena);

%% Channel Encoding
channel_coded_bits_Lena = Encoding_hamming(bits_Lena_with_CRC);

%% Interleaving
shuffled_bits_Lena = interleaving_bits(channel_coded_bits_Lena);

%% Modulation
symbols_Lena = 2*shuffled_bits_Lena - 1; % BPSK mapping

%% OFDM
symbols_time_Lena = OFDM(symbols_Lena);

%% Add Pilot Signal
[pilot_time_Lena, pilot_freq_Lena] = make_pilot(symbols_Lena);

%% Parallel to Serial
% Note that, the preamble signal is added only at first time
% And the pilot signal is added periodically every four OFDM blocks
tx_signal = Make_Tx(preamble, pilot_time_Lena, symbols_time_Lena);


%% RX_signal
rx_signal = tx_signal;

%% Time Synchronization
rx_sync = time_sync(rx_signal);

%% OFDM 1
demod_OFDM_blks = OFDM1(symbols_Lena, rx_sync);

%% Channel Estimation & Equalization
symbols_eq = CE_EQ(pilot_freq_Lena, demod_OFDM_blks);

%% Detection
% Symbol detection after the equalization
symbols_detect = symbol_detection(symbols_eq);

%% Demodulation
shuffled_bits = Demodulation(symbols_detect);

%% Deinterleaving
decoded_bits_coded = deinterleaving(shuffled_bits);

%% Channel Decoding
[HD_decoded_bits, HD_decoded_bits_with_CRC] = Decoding_Hamming(decoded_bits_coded);

%% Check CRC
CRC_check_res = Check_CRC(HD_decoded_bits_with_CRC);

%% Find BLER
BLER = Check_BLER(CRC_check_res);

%% Find BER
[~,HD_BER] = biterr(bits_Lena,HD_decoded_bits.');

%% TX Functions
function image_bits = imagetoBits(img)
imresize_scale = 0.5; % Image scaling factor
resized_img = imresize(img, imresize_scale); % Resize the image via bicubic interpolation
gray_img = rgb2gray(resized_img); % Color to grayscale
binarized_img = imbinarize(gray_img); % Grayscale to monochrome
image_bits = binarized_img(:); % Matrix vectorization
end

function bits_with_CRC = ADD_CRC(image_bits)
poly = [1 0 1 1];
reshaped_btis = reshape(image_bits,8,length(image_bits)/8);
reshaped_transposed_bits = reshaped_btis.';
add_zero_bits = [reshaped_transposed_bits zeros(height(reshaped_transposed_bits),3)];

for i = 1:height(add_zero_bits)
    remainder_Tx = [reshaped_transposed_bits(i,:) 0 0 0];
    for j = 1:size(reshaped_transposed_bits,2)
        if remainder_Tx(j) == 1
            remainder_Tx(j:j+length(poly)-1) = xor(remainder_Tx(j:j+length(poly)-1),poly);
        end
    end
    CRC = remainder_Tx(end-2:end);
    add_zero_bits(i,:) = [reshaped_transposed_bits(i,:) CRC];
end

b = add_zero_bits.';
bits_with_CRC = b(:);
end

function channel_coded_bits = Encoding_hamming(bits)
k = 11; % Number of valid bits in 1 coded bits
c = 4; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
P = [1 0 1 1;1 1 1 0;1 1 0 1;1 1 0 0;0 1 1 1;1 0 1 1;1 1 1 1;0 1 0 1;1 0 1 0;0 1 1 0;1 0 0 1]; % Parity generator matrix
G = [eye(k),P]; % code generator matrix

reshaped_bits = reshape(bits,[k,length(bits)/k]);
haming_coded_bits = mod(transpose(reshaped_bits)*G,2);

transpose_haming_coded_bits = haming_coded_bits.';
channel_coded_bits = transpose_haming_coded_bits(:);
end

function shuffled_bits = interleaving_bits(channel_coded_bits)
rng('default');
% 무작위 인덱스를 생성
numElements = numel(channel_coded_bits);
randomIndex = randperm(numElements);

% intrlv 함수를 이용하여 순서 바꾸기
shuffled_bits = intrlv(channel_coded_bits, randomIndex);
end

function symbols_time = OFDM(symbols)
% Modulation & Parameter Setting
M = length(symbols); % Number of bits 
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = M/(N/4); % Number of valid OFDM blocks
N_blk = cn + cn/4; % Number of OFDM blocks including pilot signal

% Serial to Parallel
symbols_freq={};
for i = 1:cn
    symbols_freq{end+1} = [zeros(N/4,1);0;symbols(N/4*(i-1)+1:N/4*i)]; % 64개만 사용
    symbols_freq{end} = [symbols_freq{end}; flip(symbols_freq{end}(2:end-1))];
end

% Inverse Discrete Fourier Transform (IDFT)
symbols_time={};
for i = 1:length(symbols_freq)
    symbols_time{end+1} = ifft(symbols_freq{i},N) * sqrt(N);
end

% Insert Cyclic Prefix
for i = 1:length(symbols_time)
    symbols_time{i}=[symbols_time{i}(end-N_cp+1:end); symbols_time{i}];
end

end

function tx_signal = Make_Tx(preamble,pilot_time,symbols_time)
tx_signal = [preamble; pilot_time];
for i = 1:length(symbols_time)
    tx_signal = [tx_signal; symbols_time{i}];
    if rem(i,4) == 0 && i ~= length(symbols_time)
        tx_signal = [tx_signal;pilot_time];
    end
end

end

function [pilot_time, pilot_freq] = make_pilot(symbols)

M = length(symbols); % Number of bits 
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = M/(N/4); % Number of valid OFDM blocks
N_blk = cn + cn/4; % Number of OFDM blocks including pilot signal


rng('default')
pilot_half = [zeros(N/4,1);1; 2*randi([0,1],N/4,1)-1];
pilot_freq = [pilot_half; flip(pilot_half(2:end-1))];
pilot_time = ifft(pilot_freq)*sqrt(N);
pilot_time =[pilot_time(end-N_cp+1:end); pilot_time];
end

%% RX Functions
function rx_sync = time_sync(rx_signal)
omega = 10;
mu =0.1;
Tp = 100;
tp = (1:Tp).';
preamble = cos(omega*tp+mu*tp.^2/2);

[xC, lags] = xcorr(rx_signal, preamble);
[~,idx] = max(xC);
start_pt = lags(idx);

rx_sync = rx_signal(start_pt+Tp+1:end);

end

function demod_OFDM_blks = OFDM1(symbols, rx_signal)

M = length(symbols); % Number of bits
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = M/(N/4); % Number of valid OFDM blocks
N_blk = cn + cn/4; % Number of OFDM blocks including pilot signal

% Serial to Parallel
OFDM_blks={};
for i = 1:N_blk
    OFDM_blks{end+1} = rx_signal(N_cp+1:N+N_cp);
    rx_signal = rx_signal(N_cp+N+1:end);
end

% Discrete Fourier Transform (DFT)
demod_OFDM_blks = {};
for i = 1:length(OFDM_blks)
    demod_OFDM_blks{end+1} = fft(OFDM_blks{i})/sqrt(N); % 256 point DFT
end

end

function symbols_eq = CE_EQ(pilot_freq,demod_OFDM_blks)

symbols_eq = {};
for i = 1:length(demod_OFDM_blks)
    if rem(i,5) == 1
        channel = demod_OFDM_blks{i} ./ pilot_freq;
    else
        symbols_eq{end+1} = demod_OFDM_blks{i} ./ channel;
    end
end

end

function symbols_detect = symbol_detection(symbols_eq)

symbols_detect = {};
for i = 1:length(symbols_eq)
    symbols_detect{end+1} = sign(real(symbols_eq{i}));
end

end

function shuffled_bits = Demodulation(symbols_detect)

N = 256; % Number of subcarriers

symbols_est = [];
for i = 1:length(symbols_detect)
    symbols_est = [symbols_est; symbols_detect{i}(N/4+2:2*N/4+1)];
end
symbols_est(symbols_est==0)=1;

shuffled_bits = (symbols_est+1) / 2 ;
end

function decoded_bits_coded = deinterleaving(shuffled_bits)
rng('default');
% 무작위 인덱스를 생성
numElements = numel(shuffled_bits);
randomIndex = randperm(numElements);

% intrlv 함수를 이용하여 순서 바꾸기
decoded_bits_coded = deintrlv(shuffled_bits, randomIndex);
end

function [HD_decoded_bits, HD_decoded_bits_with_CRC] = Decoding_Hamming(decoded_bits_coded)

k = 11; % Number of valid bits in 1 coded bits
c = 4; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
P = [1 0 1 1;1 1 1 0;1 1 0 1;1 1 0 0;0 1 1 1;1 0 1 1;1 1 1 1;0 1 0 1;1 0 1 0;0 1 1 0;1 0 0 1]; % Parity generator matrix
G = [eye(k),P]; % code generator matrix

H =[P' eye(c)]; % Decoding matrix
reshaped_demodulated_bits = reshape(decoded_bits_coded, [n,length(decoded_bits_coded)/n]).';
syndrome_matrix = mod(reshaped_demodulated_bits*H.',2);

HD_decoded_bits_with_CRC = [];
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
        HD_decoded_bits_with_CRC = [HD_decoded_bits_with_CRC reshaped_demodulated_bits(ii,1:11)]; %error가 발생하지 않은 경우, 받은 비트 그대로 사용
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:8)];
    else
        reshaped_demodulated_bits(ii,find) = mod((reshaped_demodulated_bits(ii,find)+1),2); % 에러가 발생한 곳 correction
        HD_decoded_bits_with_CRC = [HD_decoded_bits_with_CRC reshaped_demodulated_bits(ii,1:11)];
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:8)];
    end
end

end

function CRC_check_res = Check_CRC(HD_decoded_bits)
poly = [1 0 1 1];
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

end

function BLER = Check_BLER(CRC_check_res)
flag = 0;
for i = 1:height(CRC_check_res)
    if CRC_check_res(i,:) ~= [0 0 0]
        flag = flag + 1;
    end
end
BLER = flag/height(CRC_check_res);
end

