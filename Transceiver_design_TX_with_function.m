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
pilot_time_Lena = make_pilot(symbols_Lena);

%% Parallel to Serial
% Note that, the preamble signal is added only at first time
% And the pilot signal is added periodically every four OFDM blocks

tx_signal = Make_Tx(preamble, pilot_time_Lena, symbols_time_Lena);

%% Transmission
fs= 10000;
disp('Transmit Signal...')
sound(tx_signal,fs) % Sampling rate : 10,000Hz

signal_duration = 20;
pause(signal_duration)


%% Record the Sound for ACK
% First, create and audioDeviceReader system object
delaytime = 5;
pause(delaytime)

devicereader = audioDeviceReader(10000);
setup(devicereader);

disp('Recording. . .')
tic; % set the timer
rx_signal = [];

while toc < 10
    acquiredAudio = devicereader();
    rx_signal = [rx_signal; acquiredAudio];
end

disp('Recording Completed')


%% Tx Functions
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

function pilot_time = make_pilot(symbols)

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

