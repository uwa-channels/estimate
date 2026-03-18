% Demonstration of the channel estimation algorithm.
%
% This script demonstrates:
%   1. Generating a synthetic multipath channel with Doppler.
%   2. Estimating the channel using the est() function.
%   3. Saving the channel in the library format (params, meta).
%   4. Plotting the results using plots().
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
%
% License: MIT
%
% Revision history:
%   - Jun.  1, 2025: Initial release.
%   - Mar.  9, 2026: Updated to use params/meta struct.
%

%% Clean
clear;
clc;
close all;

%% Parameters
Fs = 48e3;            % Sampling rate
fc = 10e3;            % Center frequency
R = 8e3;              % Symbol rate
Ns = Fs / R;          % Samples per symbol
M = 4;                % Number of receivers

a = -1 / 1500;        % a = v / c
Tmp = 10e-3;          % Multipath spread
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

%% Estimation input
data.Fs = Fs;
data.R = R;
data.fc = fc;
data.d = d_channel;
data.u = u_all_channel;
data.r = received_channel;
data.element_spacing = 0.075;
data.vertical = true;

%% Channel estimation parameters
param = struct;
param.nsd = 2;                  % Samples per symbol in delay domain
param.nslr = 2;                 % Delay tracking rate: Tr = Ts / nslr
param.nst = 1 / 100;            % Samples per symbol in time domain
param.M = size(data.r, 2);
param.N = 1;                    % Number of preambles to estimate
param.K_1 = 30;                 % Anti-causal symbols (excluding 0 delay)
param.K_2 = 100;                % Causal symbols (including 0 delay)
param.delta = 0;                % Synchronization fine adjustment [samples]
param.Fs = data.Fs;
param.R = data.R;
param.fc = data.fc;
param.element_spacing = data.element_spacing;
param.vertical = data.vertical;
param.sync_length = 255;
param.limit = -30;

param.optim = 1;                % 1: LMS, 2: RLS, 3: SFTF
if param.optim == 1
  param.mu = 0.3 / (param.K_1 + param.K_2) / param.nsd;
elseif param.optim == 2
  param.lambda = 0.995;
  param.regularization = 1e-3;
elseif param.optim == 3
  param.lambda = 0.999;
  param.regularization = 1e-3;
  param.conversion = 1;
else
  error("Wrong optimizer option.");
end

param.Kf_1 = 0.001;
param.Kf_2 = param.Kf_1 / 10;
param.delay_tracking = true;

%% Estimate
tic
[h_hat, theta_hat] = est(data, param);
toc

%% Build output structs
params = struct;
params.fs_delay = param.R * param.nsd;
params.fs_time = param.R * param.nst;
params.fc = param.fc;

meta = struct;
meta.description = sprintf('Synthetic; M=%d; v/c=%.4f; Tmp=%.0fms; fc=%.0fkHz', ...
  M, a, Tmp*1e3, fc/1e3);
meta.nsd = param.nsd;
meta.nst = param.nst;
meta.Fs = param.Fs;
meta.K_1 = param.K_1;
meta.K_2 = param.K_2;
meta.fc = param.fc;
meta.element_spacing = param.element_spacing;
meta.vertical = param.vertical;
meta.delay_tracking = param.delay_tracking;
meta.limit = param.limit;
meta.optim = param.optim;
if meta.optim == 1
  meta.mu = param.mu;
elseif meta.optim == 2
  meta.lambda = param.lambda;
  meta.regularization = param.regularization;
elseif meta.optim == 3
  meta.lambda = param.lambda;
  meta.regularization = param.regularization;
  meta.conversion = param.conversion;
end
meta.Kf_1 = param.Kf_1;
meta.Kf_2 = param.Kf_2;
if param.delay_tracking
  meta.nslr = param.nslr;
end
meta.codename = 'synthetic_example';

%% Save
save_path = 'synthetic_example';
version = 1.0;
if param.delay_tracking
  phi_hat = theta_hat;
  save(save_path + ".mat", "h_hat", "phi_hat", "params", "meta", "version", ...
    "-v7.3", "-nocompression");
else
  save(save_path + ".mat", "h_hat", "theta_hat", "params", "meta", "version", ...
    "-v7.3", "-nocompression");
end
fprintf('Saved to %s.mat\n', save_path);

%% ====================================================================
%  Plotting
% =====================================================================
%
% Option A: Plot directly after estimation (as above).
%
%   figs = plots(h_hat, theta_hat, params, meta);
%
% Option B: Plot with standalone subplot figures for a paper.
%
%   figs = plots(h_hat, theta_hat, params, meta, true);
%   exportgraphics(figs.main, save_path + ".pdf", "Resolution", 300);
%   exportgraphics(figs.ir, save_path + "_ir.pdf", "Resolution", 300);
%
% Option C: Load a previously saved .mat file and plot.
%
%   channel = load('synthetic_example.mat');
%   if isfield(channel, 'phi_hat')
%     phase = channel.phi_hat;
%   else
%     phase = channel.theta_hat;
%   end
%   figs = plots(channel.h_hat, phase, channel.params, channel.meta);
%

% --- Plot from the live estimation ---
if param.delay_tracking
  phase = theta_hat;
else
  phase = theta_hat;
end
figs = plots(h_hat, phase, params, meta, true);

% Export
names = fieldnames(figs);
for k = 1:length(names)
  set(findall(figs.(names{k}), '-property', 'FontName'), 'FontName', 'SansSerif');
  if strcmp(names{k}, 'main')
    exportgraphics(figs.(names{k}), save_path + ".pdf", "Resolution", 300);
  else
    exportgraphics(figs.(names{k}), sprintf("%s_%s.pdf", save_path, names{k}), "Resolution", 300);
  end
end
fprintf('Figures saved to %s*.pdf\n', save_path);

%% --- Example: reload and plot from saved file ---
% Uncomment to test:
%
% channel = load('synthetic_example.mat');
% if isfield(channel, 'phi_hat')
%   phase = channel.phi_hat;
% else
%   phase = channel.theta_hat;
% end
% figs = plots(channel.h_hat, phase, channel.params, channel.meta);

%% Subfunctions
function p = pwr(x)
p = mean(abs(x).^2, 1);
end

% [EOF]
