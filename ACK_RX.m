clc
clear all
close all

ACK_bits = 1;

ACK_rep = [ACK_bits ACK_bits ACK_bits];

%% Modulation & Parameter Setting
symbols = 2*ACK_rep - 1; % BPSK mapping

M = length(symbols); % Number of bits
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = 1; % Number of valid OFDM blocks
N_blk = cn + 1; % Number of OFDM blocks including pilot signal

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

