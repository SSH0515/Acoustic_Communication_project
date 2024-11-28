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

%% Channel Encoding
k = 4; % Number of valid bits in 1 coded bits
c = 3; % Numbr of parities in 1 coded bits
n = k + c; % Number of total bits in 1 coded bits
A = [1 1 1;1 1 0;1 0 1;0 1 1]; % Parity generator matrix
G = [eye(k),A]; % code generator matrix

reshaped_bits = reshape(bits,[k,length(bits)/k]);
haming_coded_bits = mod(transpose(reshaped_bits)*G,2);

transpose_haming_coded_bits = haming_coded_bits.';
channel_coded_bits = transpose_haming_coded_bits(:);

%% Interleaving
rng('default');
% 무작위 인덱스를 생성
numElements = numel(channel_coded_bits);
randomIndex = randperm(numElements);

% intrlv 함수를 이용하여 순서 바꾸기
shuffled_bits = intrlv(channel_coded_bits, randomIndex);

%% Modulation & Parameter Setting
symbols = 2*shuffled_bits - 1; % BPSK mapping

M = length(symbols); % Number of bits 
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = M/(N/4); % Number of valid OFDM blocks
N_blk = cn + cn/4; % Number of OFDM blocks including pilot signal

%% Serial to Parallel
% Arrange the symbols in frequency domain
% To make the result of IDFT as a real number, symmetry is given
% The symmetry is given by adding 256 subcarriers in one cell

symbols_freq={};
for i = 1:cn
    symbols_freq{end+1} = [zeros(N/4,1);0;symbols(N/4*(i-1)+1:N/4*i)]; % 64개만 사용
    symbols_freq{end} = [symbols_freq{end}; flip(symbols_freq{end}(2:end-1))];
end

%% Inverse Discrete Fourier Transform (IDFT)
% To change the symbols in the frequency domaion to the time domain

symbols_time={};
for i = 1:length(symbols_freq)
    symbols_time{end+1} = ifft(symbols_freq{i},N) * sqrt(N);
end

%% Insert Cyclic Prefix
% Add cyclic prefix to prevent ISI

for i = 1:length(symbols_time)
    symbols_time{i}=[symbols_time{i}(end-N_cp+1:end); symbols_time{i}];
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

%% Preamble
omega = 10;
mu =0.1;
Tp = 100;
tp = (1:Tp).';
preamble = cos(omega*tp+mu*tp.^2/2);

%% Parallel to Serial
% Note that, the preamble signal is added only at first time
% And the pilot signal is added periodically every four OFDM blocks

tx_signal = [preamble; pilot_time];
for i = 1:length(symbols_time)
    tx_signal = [tx_signal; symbols_time{i}];
    if rem(i,4) == 0 && i ~= length(symbols_time)
        tx_signal = [tx_signal;pilot_time];
    end
end

%% Transmission
fs= 10000;
disp('Transmit Signal...')
sound(tx_signal,fs) % Sampling rate : 10,000Hz

signal_duration = 18;
pause(signal_duration)

%% ACK parameter
ACK_bits = 1;

ACK_rep = [ACK_bits ACK_bits ACK_bits];

%% Modulation & Parameter Setting
symbols = 2*ACK_rep - 1; % BPSK mapping

M = length(symbols); % Number of bits
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = 1; % Number of valid OFDM blocks
N_blk = cn + 1; % Number of OFDM blocks including pilot signal


%% Record the Sound for ACK
% First, create and audioDeviceReader system object

devicereader = audioDeviceReader(10000);
setup(devicereader);

disp('Recording. . .')
tic; % set the timer
rx_signal = [];

while toc < 5
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

%% Channel Estimation & Equalization

symbols_eq = {};
channel = demod_OFDM_blks{1} ./ pilot_freq;
symbols_eq{end+1} = demod_OFDM_blks{2} ./ channel;

%% Detection
% Symbol detection after the equalization

symbols_detect = {};
for i = 1:length(symbols_eq)
    symbols_detect{end+1} = sign(real(symbols_eq{i}));
end

%% Demodulation
symbols_est = symbols_detect{1}(N/2-1 : N/2+1);

symbols_est(symbols_est==0)=1;

demodulated_bits = (symbols_est+1) / 2 ;

%% Channel decoding

if length(find(demodulated_bits==1)) > length(find(demodulated_bits==0))
    decoded_bits = 1;
else
    decoded_bits = 0;
end

decoded_bits


