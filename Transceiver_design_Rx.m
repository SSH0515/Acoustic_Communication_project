clc
clear all
close all

%% Source recording
imresize_scale = 0.5; % Image scaling factor
img = imread('Lena_color.png'); % Load image
resized_img = imresize(img, imresize_scale); % Resize the image via bicubic interpolation
gray_img = rgb2gray(resized_img); % Color to grayscale
binarized_img = imbinarize(gray_img); % Grayscale to monochrome
bits = binarized_img(:); % Matrix vectorization

%% Modulation & Parameter Setting
symbols = 2*bits - 1; % BPSK mapping

M = 7/4*length(symbols); % Number of bits
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = M/(N/4); % Number of valid OFDM blocks
N_blk = cn + cn/4; % Number of OFDM blocks including pilot signal

%% Preamble
omega = 10;
mu =0.1;
Tp = 100;
tp = (1:Tp).';
preamble = cos(omega*tp+mu*tp.^2/2);

%% Record the Sound for 10 seconds
% First, create and audioDeviceReader system object
devicereader = audioDeviceReader(10000);
setup(devicereader);

disp('Recording. . .')
tic; % set the timer
rx_signal = [];

while toc < 20
    acquiredAudio = devicereader();
    rx_signal = [rx_signal; acquiredAudio];
end

disp('Recording Completed')

%% Time Synchronization
[xC, lags] = xcorr(rx_signal, preamble);
[~,idx] = max(xC);
start_pt = lags(idx);

rx_signal = rx_signal(start_pt+Tp+1:end);

%% Serial to Parallel
% Delete the cyclic prefix and make it parallel

OFDM_blks={};
for i = 1:N_blk
    OFDM_blks{end+1} = rx_signal(N_cp+1:N+N_cp);
    rx_signal = rx_signal(N_cp+N+1:end);
end

%% Discrete Fourier Transform (DFT)
% To change the symbols in the time domain to the frequency domain

demod_OFDM_blks = {};
for i = 1:length(OFDM_blks)
    demod_OFDM_blks{end+1} = fft(OFDM_blks{i})/sqrt(N); % 256 point DFT
end

%% Add Pilot Signal
% Pilot signal for channel estimation
% Firstly, generate it in the frequency domain
% Then, change it to the time domain by using IDFT
% The cyclic prefix is also added to pilot signal

rng('default')
pilot_half = [zeros(N/4,1);1; 2*randi([0,1],N/4,1)-1];
pilot_freq = [pilot_half; flip(pilot_half(2:end-1))];
pilot_time = ifft(pilot_freq)*sqrt(N);
pilot_time =[pilot_time(end-N_cp+1:end); pilot_time];


%% Channel Estimation & Equalization

symbols_eq = {};
for i = 1:length(demod_OFDM_blks)
    if rem(i,5) == 1
        channel = demod_OFDM_blks{i} ./ pilot_freq;
    else
        symbols_eq{end+1} = demod_OFDM_blks{i} ./ channel;
    end
end

%% Detection
% Symbol detection after the equalization

symbols_detect = {};
for i = 1:length(symbols_eq)
    symbols_detect{end+1} = sign(real(symbols_eq{i}));
end

%% Demodulation
symbols_est = [];
for i = 1:length(symbols_detect)
    symbols_est = [symbols_est; symbols_detect{i}(N/4+2:2*N/4+1)];
end
symbols_est(symbols_est==0)=1;

shuffled_bits = (symbols_est+1) / 2 ;

%% Deinterleaving
rng('default');
% 무작위 인덱스를 생성
numElements = numel(shuffled_bits);
randomIndex = randperm(numElements);

% intrlv 함수를 이용하여 순서 바꾸기
decoded_bits_coded = deintrlv(shuffled_bits, randomIndex);

%% Variables
k = 4; % Number of valid bits in 1 coded bits
c = 3; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
A = [1 1 1;1 1 0;1 0 1;0 1 1]; % Parity generator matrix
G = [eye(k),A]; % code generator matrix

bit_matrix = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1; ...
        0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1; ...
        1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1; ...
        1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1];

code_matrix = mod(bit_matrix*G,2); % code matrix

%% Channel Decoding
demodulated_bits = decoded_bits_coded; % BPSK Demodulation

H =[A' eye(c)]; % Decoding matrix
reshaped_demodulated_bits = reshape(demodulated_bits, [n,length(demodulated_bits)/n]).';
syndrome_matrix = mod(reshaped_demodulated_bits*H.',2);

error = 0; % Block error rate
for hh = 1 : height(reshaped_demodulated_bits)
    hamming_distance = zeros(height(code_matrix),1);
    for num = 1:height(code_matrix)  
        hamming_distance(num) = sum(reshaped_demodulated_bits(hh,:) ~= code_matrix(num,:));
    end
    if min(hamming_distance) >= 2
        error = error+1;
    end
end
BLER = error/height(reshaped_demodulated_bits); % Calculation BLER

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
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:4)]; %error가 발생하지 않은 경우, 받은 비트 그대로 사용
    else
        reshaped_demodulated_bits(ii,find) = mod((reshaped_demodulated_bits(ii,find)+1),2); % 에러가 발생한 곳 correction
        HD_decoded_bits = [HD_decoded_bits reshaped_demodulated_bits(ii,1:4)];
    end
end

[~,HD_BER] = biterr(bits,HD_decoded_bits);

%% Source Decoding & Show Image
estimated_img = reshape(HD_decoded_bits , [sqrt(length(HD_decoded_bits)),sqrt(length(HD_decoded_bits))]);
resized_estimated_img = imresize(estimated_img, 1/imresize_scale);
imshow(resized_estimated_img)