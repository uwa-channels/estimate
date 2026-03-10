function figs = plots(h_hat, phase, params, meta, subplots)
% PLOTS Visualize the channel impulse responses.
%
% FIGS = PLOTS(H_HAT, PHASE, PARAMS, META) creates the 2x3 summary
% figure and returns it in FIGS.main.
%
% FIGS = PLOTS(H_HAT, PHASE, PARAMS, META, TRUE) also creates
% standalone subplot figures.  FIGS is a struct with fields:
%   main         - 2x3 summary figure
%   ir           - Time-varying impulse response
%   doppler      - Doppler-delay
%   delay_angle  - Delay-angle (if vertical array with M >= 13)
%   mip          - Multipath intensity profile
%   pll          - Phase-locked loop
%
% Saving is left to the caller, e.g.:
%   figs = plots(h_hat, theta_hat, params, meta, true);
%   exportgraphics(figs.main, "for_paul/yellow_2.pdf");
%   exportgraphics(figs.ir, "for_paul/yellow_2_ir.pdf");
%
% Inputs:
%    h_hat          - Channel impulse response, size (K, M, T).
%    phase          - Phase trajectory (theta_hat or phi_hat), size (M, N).
%    params         - Struct with fields: fs_delay, fs_time, fc.
%    meta           - Struct with estimation metadata (K_1, K_2, nsd,
%                     nst, fc, element_spacing, vertical, optim, mu or
%                     lambda, Kf_1, Kf_2, delay_tracking, codename,
%                     and optionally nslr).
%    subplots       - Logical, create standalone figures (default false).
%
% Example:
%    load('yellow_2.mat');
%    figs = plots(h_hat, theta_hat, params, meta);
%
% License: MIT
%
% Authors: Zhengnan Li, Diego A. Cuji
% Emails : uwa-channels@ofdm.link, cujidutan.d@northeastern.edu
%
% Revision history:
%   - Jul. 09, 2024: Initial release; add delay-angle plot.
%   - Jul. 11, 2024: Add capabilities of resampling the time domain.
%   - Mar. 09, 2026: Refactored to use params/meta struct. Returns
%                     struct of figure handles.

%% Font settings
set(groot, 'defaultAxesFontName', 'SansSerif');
set(groot, 'defaultTextFontName', 'SansSerif');

%% Default arguments
if nargin < 5, subplots = false; end

%% Derive parameters from params + meta
K_1 = meta.K_1;
K_2 = meta.K_2;
M = size(h_hat, 2);
ns = meta.nsd;
fc = params.fc;
d = meta.element_spacing;
R = params.fs_delay / ns;
fs = params.fs_delay;
delay_tracking = meta.delay_tracking;

% Plotting lower limit (dB)
if isfield(meta, 'limit')
  limit = meta.limit;
else
  limit = -30;
end

% Channel name for title
if isfield(meta, 'codename')
  name = meta.codename;
else
  name = 'channel';
end

%% Colors and figure settings
% colors = flipud([0.9059    0.3843    0.3294;
%   0.9373    0.5412    0.2784;
%   0.9686    0.6667    0.3451;
%   1.0000    0.8157    0.4353;
%   1.0000    0.9020    0.7176;
%   0.6667    0.8627    0.8784;
%   0.4471    0.7373    0.8353;
%   0.3216    0.5608    0.6784;
%   0.2157    0.4039    0.5843;
%   0.1176    0.2745    0.4314]);
colors = [1.0000    1.0000    1.0000;
  0.9176    0.9490    0.9804;
  0.7765    0.8588    0.9373;
  0.6196    0.7647    0.8863;
  0.4549    0.6627    0.8118;
  0.2863    0.5451    0.7451;
  0.1569    0.4431    0.6745;
  0.0706    0.3451    0.5804;
  0.0392    0.2627    0.4863;
  0.0314    0.2039    0.4039;
  0.0196    0.1490    0.3137;
  0.0118    0.0941    0.2157;
  0.0039    0.0392    0.1098];
cmap = interp1(colors, linspace(1, size(colors,1), 64));
linecolors = interp1(colors, linspace(1, size(colors,1), M));
figs.main = figure('position', [0, 0, 1200, 900]);
colormap(cmap);

%% Some helper variables
K = K_1 + K_2;
if size(h_hat, 2) > 1
  h_hat_first = squeeze(h_hat(:, floor(size(h_hat, 2)/2), :)).';
else
  h_hat_first = squeeze(h_hat(:, 1, :)).';
end
h_hat_first_channel = h_hat_first.';
h_hat_first_channel = h_hat_first_channel ./ max(abs(h_hat_first_channel), [], "all");

%% Axis variables
if R > 1000
  xaxis = (-K_1*ns:K_2*ns-1) / fs * 1e3;
  xunit = 'Delay [ms]';
else
  xaxis = (-K_1*ns:K_2*ns-1) / fs;
  xunit = 'Delay [s]';
end
yaxis = (1:size(h_hat_first_channel,2)) / params.fs_time;

%% Subplot 1, time varying impulse response
set(0, 'currentfigure', figs.main);
ax = subplot(231);
imagesc(xaxis, yaxis, 20*log10(abs(h_hat_first_channel)).', [limit, 0]);
xlabel(xunit), ylabel('Time [s]');
c = colorbar(ax, "Location", "Northoutside");
ps = get(ax, "Position"); pc = get(c, "Position"); pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; ps(4)=ps(4)*0.8;set(c, "Position", pc); set(ax, "Position", ps); title(c, "dB"); lb = get(c, "Title"); set(lb, "Units", "Normalized", "Position", [1.15, 0]);
if delay_tracking
  text(ax, 0.55, 1.09, 0.0, "$\underline{\hat{\mathbf{h}}}[n]$", "Interpreter", "latex", "Units", "Normalized")
else
  text(ax, 0.55, 1.09, 0.0, "$\hat{\mathbf{h}}[n]$", "Interpreter", "latex", "Units", "Normalized")
end

if subplots
  figs.ir = figure('position', [0, 0, 400, 400]); ax = gca;
  imagesc(xaxis, yaxis, 20*log10(abs(h_hat_first_channel)).', [limit, 0]);
  xlabel(xunit), ylabel('Time [s]');
  c = colorbar(ax, "Location", "Northoutside");
  ps = get(ax, "Position"); ps(4)=ps(4)*0.9; set(ax, "Position", ps);
  pc = get(c, "Position"); pc(2)=pc(2)*1.03; pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; set(c, "Position", pc);
  title(c, "dB"); lb_sub = get(c, "Title"); set(lb_sub, "Units", "Normalized", "Position", [1.15, 0]);
  if delay_tracking
    text(ax, 0.65, 1.09, 0.0, "$\underline{\hat{\mathbf{h}}}[n]$", "Interpreter", "latex", "Units", "Normalized")
  else
    text(ax, 0.65, 1.09, 0.0, "$\hat{\mathbf{h}}[n]$", "Interpreter", "latex", "Units", "Normalized")
  end
  colormap(cmap);
end

%% Subplot 2, Doppler-delay (scattering function)
% Blackman-Tukey estimator: S(τ,ν) = FFT{ R_h(τ, Δt) }
set(0, 'currentfigure', figs.main);
ax = subplot(232);
N_time = size(h_hat_first_channel, 2);
nfft = 2^(nextpow2(2*N_time - 1) + 1);
doppler_spectrum = zeros(size(h_hat_first_channel, 1), nfft);
for ii = 1:size(h_hat_first_channel, 1)
  rr = xcorr(h_hat_first_channel(ii, :));
  doppler_spectrum(ii, :) = abs(fftshift(fft(rr, nfft)));
end
delay_doppler = doppler_spectrum / max(doppler_spectrum, [], 'all');
yaxis_doppler = linspace(-params.fs_time/2, params.fs_time/2 - params.fs_time/nfft, nfft);
imagesc(xaxis, yaxis_doppler, 10*log10(delay_doppler).', [limit, 0]);
xlabel(xunit); ylabel('Doppler [Hz]');
ylim([-10, 10])
c = colorbar(ax, "Location", "Northoutside");
ps = get(ax, "Position"); pc = get(c, "Position"); pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; ps(4)=ps(4)*0.8; set(c, "Position", pc); set(ax, "Position", ps); title(c, "dB"); lb = get(c, "Title"); set(lb, "Units", "Normalized", "Position", [1.15, 0]);
text(ax, 0.55, 1.09, 0.0, "$\left|\hat{S}_{\hat{\underline{h}}}(\tau, \nu)\right|$", "Interpreter", "latex", "Units", "Normalized")

if subplots
  figs.doppler = figure('position', [0, 0, 400, 400]); ax = gca;
  imagesc(xaxis, yaxis_doppler, 10*log10(delay_doppler).', [limit, 0]);
  xlabel(xunit); ylabel('Doppler [Hz]');
  ylim([-2, 2])
  c = colorbar(ax, "Location", "Northoutside");
  ps = get(ax, "Position"); ps(4)=ps(4)*0.9; set(ax, "Position", ps);
  pc = get(c, "Position"); pc(2)=pc(2)*1.0; pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; set(c, "Position", pc);
  title(c, "dB"); lb_sub = get(c, "Title"); set(lb_sub, "Units", "Normalized", "Position", [1.15, 0]);
  text(ax, 0.65, 1.09, 0.0, "$\left|\hat{S}_{\hat{\underline{h}}}(\tau, \nu)\right|$", "Interpreter", "latex", "Units", "Normalized")
  colormap(cmap);
end

%% Subplot 3, Array
set(0, 'currentfigure', figs.main);
nfft = 2^(nextpow2(K*ns));
if size(h_hat, 3) > 64
  obs = 64;
else
  obs = size(h_hat, 3);
end
nangle = 256;
c_sound = 1500;

fk = fc + linspace(-fs/2, fs/2-fs/nfft, nfft);
delays = linspace(-K_1*ns/fs, (K_2*ns-1)/fs, K*4);
theta = linspace(-pi/4, pi/4, nangle);
exponents = exp(2j * pi * delays.' * (fk-fc));

% Step 1: formulate the broadband beamformer
s_prime = zeros(nfft, M, nangle);
for m = 1 : M
  s_prime(:, m, :) = exp(2j * pi * fk.' * d / c_sound * (m-1) * sin(theta));
end

% Step 2: convert to frequency domain, apply the beamformer, and convert back to delay domain
h_hat_bf = zeros(length(delays), nangle, obs);
for o = 1 : obs
  h_obs = h_hat(:, :, end-obs+o);
  H_hat = fft(circshift([h_obs; ...
    zeros(nfft-size(h_obs,1), M)], -K_1*ns));
  H_hat_bf = squeeze(sum(H_hat .* s_prime, 2));
  h_hat_bf(:, :, o) = exponents * H_hat_bf;
end
zaxis = sum(abs(h_hat_bf).^2, 3);
zaxis = zaxis ./ max(zaxis, [], "all");
zaxis = 10 * log10(zaxis);

if meta.vertical && M >= 13
  set(0, 'currentfigure', figs.main);
  ax = subplot(233);
  imagesc(delays * 1e3, theta/pi*180, zaxis.', [limit, 0]), colormap(cmap), axis("xy")
  xlabel(xunit); ylabel("Angle [degree]")
  c = colorbar(ax, "Location", "Northoutside");
  ps = get(ax, "Position"); pc = get(c, "Position"); pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; ps(4)=ps(4)*0.8;set(c, "Position", pc); set(ax, "Position", ps); title(c, "dB"); lb = get(c, "Title"); set(lb, "Units", "Normalized", "Position", [1.15, 0]);
  text(ax, 0.55, 1.09, 0.0, "$10 \log_{10} \left( \hat{P}(\tau, \theta) \right)$", "Interpreter", "latex", "Units", "Normalized")

  if subplots
    figs.delay_angle = figure('position', [0, 0, 400, 400]); ax = gca;
    imagesc(delays * 1e3, theta/pi*180, zaxis.', [limit, 0]), colormap(cmap), axis("xy");
    xlabel(xunit); ylabel("Angle [degree]")
    c = colorbar(ax, "Location", "Northoutside");
    ps = get(ax, "Position"); ps(4)=ps(4)*0.9; set(ax, "Position", ps);
    pc = get(c, "Position"); pc(2)=pc(2)*1.03; pc(3)=pc(3)*0.4; pc(4)=pc(4)*0.8; set(c, "Position", pc);
    title(c, "dB"); lb_sub = get(c, "Title"); set(lb_sub, "Units", "Normalized", "Position", [1.15, 0]);
    text(ax, 0.65, 1.09, 0.0, "$10 \log_{10} \left( \hat{P}(\tau, \theta) \right)$", "Interpreter", "latex", "Units", "Normalized")
    colormap(cmap);
  end
else
  set(0, 'currentfigure', figs.main);
  ax = subplot(233); hold on; view(120, 60);
  for m = 1 : M
    plot3(ax, m * ones(K*ns, 1), xaxis, abs(h_hat(:, m, end)), "Color", linecolors(m, :));
  end
  xlabel("Receiver index"); ylabel(xunit); zlabel("$\left| \hat{h}(\tau) \right|$", "Interpreter", "latex")
end

%% Subplot 4, delay profile
set(0, 'currentfigure', figs.main);
ax = subplot(234); hold on; box("on"); axis("square"); grid("on");
mip = zeros(size(h_hat, 1), M);
for m = 1:M
  mip(:, m) = mean(abs(squeeze(h_hat(:, m, :))).^2, 2);
end
mip = mip ./ max(abs(mip), [], "all");
for m = 1:M
  plot(ax, xaxis, 10*log10(abs(mip(:, m))), 'Color', linecolors(m, :));
end
xlim([min(xaxis), max(xaxis)]);
xlabel(xunit); ylabel('Multipath intensity profile [dB]');

if subplots
  figs.mip = figure('position', [0, 0, 400, 400]); ax = gca;
  hold on; box("on"); grid("on");
  for m = 1:M
    plot(ax, xaxis, 10*log10(abs(mip(:, m))), 'Color', linecolors(m, :));
  end
  xlim([min(xaxis), max(xaxis)]);
  xlabel(xunit); ylabel('Multipath intensity profile [dB]');
end

%% Subplot 5, PLL
set(0, 'currentfigure', figs.main);
ax = subplot(235); hold on; box("on"); axis("square"); grid("on");
xaxis_phase = (1:size(phase, 2)) / params.fs_delay;
for m = 1:M
  plot(ax, xaxis_phase, phase(m, :), "Color", linecolors(m, :));
end
xlabel('Time [s]'); ylabel('Phase [rad]');
xlim([min(xaxis_phase), max(xaxis_phase)]);

if subplots
  figs.pll = figure('position', [0, 0, 400, 400]); ax = gca; hold on;
  box("on"); axis("square"); grid("on");
  for m = 1:M
    plot(ax, xaxis_phase, phase(m, :), "Color", linecolors(m, :));
  end
  xlabel('Time [s]'); ylabel('Phase [rad]');
  xlim([min(xaxis_phase), max(xaxis_phase)]);
end

%% Subplot 6, params
set(0, 'currentfigure', figs.main);
ax = subplot(236); axis("off");

% Signal parameters
if R > 1000
  data_params = { ...
    '\textbf{Signal parameters}', ...
    sprintf('$F_s = %.2f$ ksa/s, $R = %.2f$ ksym/s, $f_c = %.2f$ kHz', ...
      fs/1e3, R/1e3, fc/1e3), ...
    sprintf('$M = %d$', M)};
else
  data_params = { ...
    '\textbf{Signal parameters}', ...
    sprintf('$F_s = %.2f$ sa/s, $R = %.2f$ sym/s, $f_c = %.2f$ Hz', ...
      fs, R, fc), ...
    sprintf('$M = %d$', M)};
end

% Estimation parameters
if meta.optim == 1
  step_str = sprintf('$\\mu = %s/(K_1{+}K_2)/n_{sd}$', ...
    ltxnum(meta.mu * (meta.K_1 + meta.K_2) * meta.nsd));
else
  step_str = sprintf('$\\lambda = %s$', ltxnum(meta.lambda));
end

if delay_tracking
  time_str = sprintf('$K_1 = %d$, $K_2 = %d$, $T_s = T/%d$, $T_c = T/%d$, $T_r = T/%s$', ...
    K_1, K_2, meta.nsd, round(1/meta.nst), ltxnum(meta.nsd * meta.nslr));
else
  time_str = sprintf('$K_1 = %d$, $K_2 = %d$, $T_s = T/%d$, $T_c = T/%d$', ...
    K_1, K_2, meta.nsd, round(1/meta.nst));
end

est_params = { ...
  '\textbf{Estimation parameters}', ...
  time_str, ...
  sprintf('%s, $K_{f_1} = %s$, $K_{f_2} = %s$', ...
    step_str, ltxnum(meta.Kf_1), ltxnum(meta.Kf_2))};

text(ax, 0, 0.85, data_params, "Units", "normalized", ...
  "Interpreter", "latex", "VerticalAlignment", "top", "FontSize", 10)
text(ax, 0, 0.50, est_params, "Units", "normalized", ...
  "Interpreter", "latex", "VerticalAlignment", "top", "FontSize", 10)

%%
sgtitle(name, "Interpreter", "none");

end


function s = ltxnum(x)
% LTXNUM Format a number as a LaTeX string.
%   0         => '0'
%   1, -3     => '1', '-3'
%   0.5       => '0.5'
%   1e-4      => '10^{-4}'
%   2.5e-3    => '2.5 \times 10^{-3}'
if x == 0
  s = '0';
  return
end
% Simple numbers: use plain decimal notation
if abs(x) >= 0.01 && abs(x) < 1e4
  s = num2str(x, '%.4g');
  return
end
e = floor(log10(abs(x)));
m = x / 10^e;
if abs(m - round(m)) < 1e-10
  m = round(m);
end
if m == 1
  s = sprintf('10^{%d}', e);
elseif m == -1
  s = sprintf('-10^{%d}', e);
else
  s = sprintf('%s \\times 10^{%d}', num2str(m, '%.4g'), e);
end
end
