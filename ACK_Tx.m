clc
clear all
close all

%% Channel Encoding
channel_coded_bits =zeros(3,1);
ACK =1;
if ACK == 1
    channel_coded_bits =[1;1;1];
elseif ACK == 0
    channel_coded_bits = [0;0;0];
end

%% Modulation & Parameter Setting
symbols = 2*channel_coded_bits - 1; % BPSK mapping

M = length(symbols); % Number of bits
N = 256; % Number of subcarriers
N_cp = 32; %Length of cyclic prefix
cn = 1; % Number of valid OFDM blocks
N_blk = cn + 1; % Number of OFDM blocks including pilot signal

%% Serial to Parallel
% Arrange the symbols in frequency domain
% To make the result of IDFT as a real number, symmetry is given
% The symmetry is given by adding 256 subcarriers in one cell

symbols_freq={};

symbols_freq{end+1} = [zeros(126,1);symbols(1:3,1)]; % ACK
symbols_freq{end} = [symbols_freq{end}; flip(symbols_freq{end}(2:end-1))];


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

