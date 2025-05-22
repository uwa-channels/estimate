% Demonstration of the channel estimation algorithm.
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
%
% License: MIT
%
% Revision history:
%   - Jun. 1, 2025: initial release.
%

%% Clean
clear;
clc;
close all;

%% Parameters
Fs = 48e3; % Sampling rate
fc = 10e3; % Center frequency
R = 8e3; % Symbol rate
Ns = Fs / R; % Samples per symbol
M = 4; % Number of receivers

a = -1 / 1500; % a = v / c
Tmp = 10e-3; % Multipath spread
path_delay = [0, randsample(1:Tmp*1e3-1, 8)] / 1e3;
path_gain = exp(-path_delay*1.5./Tmp) .* randsample([-1, 1], length(path_delay), true);

%% Prepare signal for channel estimation
d_channel = randsample([-1, +1], 2^17-1, true).';
u_all_channel = resample(d_channel, Ns, 1);
s_all_channel = real(u_all_channel.*exp(1j*2*pi*fc*(0:length(u_all_channel) - 1).'/Fs));
s_all_channel = [zeros(round(Fs/16), 1); s_all_channel; zeros(round(Fs/16), 1);];

%% Create multipath signal
received_channel = zeros(length(s_all_channel), M);
for m = 1:M
    for p = 1:length(path_gain)
        received_channel(:, m) = received_channel(:, m) + path_gain(p) * ...
            circshift(s_all_channel, round(path_delay(p)*Fs));
    end
end

[p, q] = rat(1-a);
received_channel = resample(received_channel, p, q, "Dimension", 1);
received_channel = received_channel ./ (sqrt(pwr(received_channel)));
received_channel = received_channel + 0.1 * randn(size(received_channel));

%% Estimate the channel
data.Fs = Fs;
data.R = R;
data.fc = fc;
data.d = d_channel;
data.u = u_all_channel;
data.r = received_channel;

%% Channel estimation parameters
params.nsd = 2; % Number of samples per symbol in delay domain
params.nslr = 2; % Tr = Ts/nslr = 1/(R*nsd*nslr) to perform the LR operation
params.nst = 1 / 100; % Number of samples per symbol in time domain
params.M = size(data.r, 2); % Number of receivers

params.N = 1; % Number of preambles to estimate
params.K_1 = 30; % Number of symbols to the left (anti-causal, excluding 0 delay)
params.K_2 = 100; % Number of symbols to the right (causal, including 0 delay)
params.delta = 0; % Synchronization fine adjustment point, in samples (Fs)
params.Fs = data.Fs; % Sample rate of the original signal
params.R = data.R; % Symbol rate in the original signal
params.fc = data.fc; % Center frequency of the original signal
params.sync_length = 255;

params.optim = 3; % 1 for LMS, 2 for RLS and 3 for SFTF
if params.optim == 1
    params.mu = 0.3 / (params.K_1 + params.K_2) / params.nsd; % LMS step size
elseif params.optim == 2
    params.lambda = 0.995; % RLS forgetting parameter
    params.regularization = 1e-3; % RLS regularization factor
elseif params.optim == 3
    params.lambda = 0.999; % SFTF forgetting parameter
    % 1e-3 less volatile; 5e-2 or more for more volatile channels
    params.regularization = 1e-3; % SFTF regularization factor;
    params.conversion = 1; % SFTF conversion factor, gamma
else
    error("Wrong optimizer option.");
end

params.Kf_1 = 0.001; % PLL coefficients 1
params.Kf_2 = params.Kf_1 / 10; % PLL coefficients 2
params.delay_tracking = true;

save_files = true;
enable_plots = true;
filename = 'test';

tic
[h_hat, theta_hat] = est(data, params);
toc

%% Save
if save_files
    params_save = params;
    params = struct;
    params.fs_delay = params_save.R * params_save.nsd;
    params.fs_time = params_save.R * params_save.nst;
    params.fc = params_save.fc;

    version = 1.0;

    file_name = sprintf("%s.mat", filename);
    if isfield(data, "f_resamp")
        f_resamp = data.f_resamp;
        save(file_name, "h_hat", "theta_hat", "params", "f_resamp", "version", ...
            "-v7.3", "-nocompression");
    else
        save(file_name, "h_hat", "theta_hat", "params", "version", ...
            "-v7.3", "-nocompression");
    end
end

%% Plots
if enable_plots
    name = sprintf("extracted/%s", filename);
    plots(h_hat, params);
end

%% Subfunctions
function p = pwr(x)
p = mean(abs(x).^2, 1);
end

function plots(h_hat, params)
delay_axis = (0:size(h_hat, 1) - 1) ./ params.fs_delay;
time_axis = (0:size(h_hat, 3) - 1) ./ params.fs_time;
figure
imagesc(delay_axis*1e3, time_axis, 20*log10(squeeze(abs(h_hat(:, 1, :))).'), [-30, 0])
xlabel('Delay [ms]')
ylabel('Time [s]')
end

% [EOF]
