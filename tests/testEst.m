classdef testEst < matlab.unittest.TestCase
  % TESTEST Unit tests for est.m and plots.m
  %
  % Each test transmits a known BPSK sequence through a known multipath
  % channel, estimates the channel with est, and checks the estimate
  % against the channel that made the signal:
  %
  %   testRecoversPaths:   h_hat equals the impulse response the channel
  %                        must produce -- every path, with its delay, gain
  %                        and carrier phase, shaped by the transmit pulse,
  %                        on every element -- and not its conjugate.
  %   testTracksDoppler:   the PLL phase slope equals -2*pi*fc*a.  With delay
  %                        tracking the taps stay put; without it they drift
  %                        by a*t.
  %   testOutputFormat:    sizes, and the duration constraint of the library
  %                        file format.
  %   testReplayRoundTrip: a channel estimated here and replayed with
  %                        uwa-channels/matlab reproduces the true channel's
  %                        output, early in the record and as late as it
  %                        allows, and does not with the sign of phi_hat
  %                        flipped.  Skipped when replay is not on the path.
  %   testRejectsInvalidInputs, testPlots.
  %
  % Only LMS (optim 1) and RLS (optim 2) are tested.  The PLL gains are
  % Kf_1 = 1e-3, Kf_2 = Kf_1/10, i.e. w_nT = 0.01 and zeta = 0.05.
  %
  % sync_length is 64 symbols.  est synchronizes with a coherent correlation
  % that does not compensate Doppler, so over 255 symbols at |a| = 1e-3 the
  % carrier rotates by about 4.8 rad and a sidelobe can beat the true peak.
  %
  % Other m-files required: est.m, plots.m; replay.m from uwa-channels/matlab
  %                         for testReplayRoundTrip.
  % Subfunctions: simulate, passband, symbol_train, default_params,
  %               expected_response, fit_response, peak_tap, compare,
  %               fit_error, plot_round_trip, lowpass_bb, pwr
  % Toolbox required: Signal Processing Toolbox.
  % MAT-files required: None.
  %
  % See also: est.m, plots.m
  %
  % Author: Zhengnan Li, Claude (Opus 5.5)
  % Email : uwa-channels@ofdm.link
  %
  % License: MIT
  %
  % Revision history:
  %   - Oct. 26, 2026: Initial release.
  %   - Oct. 27, 2026: testReplayRoundTrip also replays at the latest start
  %                     the record allows, where phi_hat has accumulated the
  %                     most drift.  It replays the symbol train rather than
  %                     the pulse-shaped probe, since h_hat holds the pulse;
  %                     compare aligns continuously; the thresholds follow;
  %                     and the round trip is logged as a figure.

  properties (Constant)
    Fs = 48e3;                          % Passband sampling rate [Hz]
    fc = 12e3;                          % Center frequency [Hz]
    R = 4e3;                            % Symbol rate [Baud]
    % Delays are whole samples at Fs, so the simulated channel is exactly the
    % one written here, and they fall a quarter, a half and three quarters of
    % a carrier cycle after a whole cycle, so every path has its own phase.
    delays = [0, 61, 166, 295] / 48e3;  % Path delays [s]
    gains = [1, -0.6, 0.4, 0.25];       % Path gains
    element_delay = 26 / 48e3;          % Extra delay per array element [s]
  end

  properties (TestParameter)
    algorithm = struct('LMS', 1, 'RLS', 2);
    doppler = struct('static', 0, 'receding', 1e-3, 'closing', -1e-3);
    drift = struct('receding', 5e-4, 'closing', -5e-4);
    tracking = struct('delay_tracking', true, 'phase_only', false);
  end

  methods (TestMethodSetup)
    function setRandomSeed(~)
      rng(1994);
    end

    function hideFigures(testCase)
      old = get(0, 'DefaultFigureVisible');
      set(0, 'DefaultFigureVisible', 'off');
      testCase.addTeardown(@() set(0, 'DefaultFigureVisible', old));
    end
  end

  methods (Test)

    function testRecoversPaths(testCase, algorithm, doppler)
      M = 2;
      [data, params] = simulate(testCase, M, doppler);
      params.optim = algorithm;
      h_hat = est(data, params);

      % Average over the last tenth of the record, where the filter has
      % converged.  Delay tracking keeps the taps in place.
      T = size(h_hat, 3);
      h = mean(h_hat(:, :, end-ceil(T/10)+1:end), 3);

      % Compare with the response the channel must produce.  The fit has
      % few free parameters: one complex gain per element, and the fraction
      % of a tap by which synchronization lands off the grid, shared by all
      % elements so that their relative delays are tested.
      [nmse, offset] = fit_response(testCase, h, params, false);
      testCase.verifyLessThan(nmse, -30, ...
        sprintf('h_hat is %.1f dB from the true impulse response.', nmse));

      % est puts the reference path K_1 symbols into the filter.
      testCase.verifyLessThanOrEqual(abs(offset), 1, ...
        sprintf('The reference path is %.2f taps from tap K_1*nsd + 1.', offset));

      % h_hat is the channel, not its conjugate.  The paths have different
      % carrier phases, so the conjugate cannot fit.
      nmse_conj = fit_response(testCase, h, params, true);
      testCase.verifyGreaterThan(nmse_conj, -10, ...
        sprintf('The conjugate of the true response fits h_hat to %.1f dB.', nmse_conj));
    end

    function testTracksDoppler(testCase, drift, tracking)
      M = 1;
      [data, params] = simulate(testCase, M, drift);
      params.optim = 2;
      params.delay_tracking = tracking;
      [h_hat, theta_hat] = est(data, params);

      fs_delay = testCase.R * params.nsd;
      fs_time = testCase.R * params.nst;

      % Phase: the received carrier runs at fc/(1+a), so the PLL phase
      % falls at 2*pi*fc*a/(1+a) rad/s.
      t = (0:size(theta_hat, 2)-1) / fs_delay;
      slope = polyfit(t, theta_hat(1, :), 1);
      expected_slope = -2*pi*testCase.fc*drift / (1+drift);
      testCase.verifyLessThan(abs(slope(1)/expected_slope - 1), 0.01, ...
        sprintf('Phase slope %.2f rad/s, expected %.2f rad/s.', ...
        slope(1), expected_slope));

      % Delay: follow the strongest path through the snapshots.  Skip the
      % first tenth, where the PLL is still pulling in.
      T = size(h_hat, 3);
      k = ceil(T/10):T;
      taps = zeros(size(k));
      for i = 1:numel(k)
        taps(i) = peak_tap(abs(h_hat(:, 1, k(i))).^2);
      end
      tap_slope = polyfit(k / fs_time, taps, 1);
      if tracking
        % The PLL phase carries the drift, so the taps must not move.
        testCase.verifyLessThan(max(taps) - min(taps), 1, ...
          sprintf('With delay tracking, the main path moved by %.2f taps.', ...
          max(taps) - min(taps)));
      else
        % The taps carry the drift: a*t seconds, i.e. a*fs_delay taps/s.
        expected = drift * fs_delay;
        testCase.verifyLessThan(abs(tap_slope(1)/expected - 1), 0.05, ...
          sprintf('The main path drifted %.3f taps/s, expected %.3f.', ...
          tap_slope(1), expected));
      end
    end

    function testOutputFormat(testCase, tracking)
      M = 3;
      [data, params] = simulate(testCase, M, 0);
      params.optim = 2;
      params.nst = 1/50;
      params.delay_tracking = tracking;
      [h_hat, theta_hat, err] = est(data, params);

      L = (params.K_1 + params.K_2) * params.nsd;
      testCase.verifySize(h_hat, [L, M, size(h_hat, 3)]);
      testCase.verifySize(theta_hat, [M, size(theta_hat, 2)]);
      testCase.verifyEqual(size(err, 2), M);
      testCase.verifyGreaterThan(size(h_hat, 3), 1);

      % The library format requires the phase and the impulse responses
      % to span the same time.
      fs_delay = testCase.R * params.nsd;
      fs_time = testCase.R * params.nst;
      testCase.verifyEqual(size(theta_hat, 2) / fs_delay, size(h_hat, 3) / fs_time, ...
        'AbsTol', 1e-12, ...
        'theta_hat and h_hat do not span the same duration.');

      testCase.verifyTrue(all(isfinite(h_hat(:))));
      testCase.verifyTrue(all(isfinite(theta_hat(:))));
    end

    function testReplayRoundTrip(testCase, doppler)
      testCase.assumeEqual(exist('replay', 'file'), 2, ...
        'replay.m from uwa-channels/matlab is not on the path.');

      M = 2;
      [data, params] = simulate(testCase, M, doppler);
      params.optim = 2;
      [h_hat, phi_hat] = est(data, params);

      % Save the estimate the way example_estimate.m does.
      channel.h_hat = h_hat;
      channel.phi_hat = phi_hat;
      channel.params.fs_delay = testCase.R * params.nsd;
      channel.params.fs_time = testCase.R * params.nst;
      channel.params.fc = testCase.fc;
      channel.version = 1.0;

      % A fresh probe, through the estimated channel and the true one.
      % h_hat maps symbols to the received signal, so it already holds the
      % transmit pulse: replay the symbols themselves, as a zero-stuffed
      % train at fs_delay carried to Fs and up to fc.  Replaying the
      % pulse-shaped probe would apply the pulse twice.
      probe = randi([0, 1], 2^13-1, 1) * 2 - 1;
      x = symbol_train(testCase, probe, params.nsd);
      truth = simulate(testCase, M, doppler, probe);
      truth = truth.r;

      % Replay early in the record and as late as it allows.  phi_hat
      % accumulates, so late in the record its absolute value is about 144
      % samples at Fs here; replay must take the drift from its starting
      % sample, or the output is shifted by all of it and loses its end
      % (receding) or its beginning (closing).
      fs_delay = channel.params.fs_delay;
      T_max = size(h_hat, 3) / channel.params.fs_time * fs_delay;
      T = ceil(numel(x) * fs_delay / testCase.Fs);
      late = floor(T_max - T - size(h_hat, 1));

      % The estimate's own noise puts the static case at about -47 dB.
      % With Doppler, est's linear interpolation in delay tracking costs
      % about 12 dB more (replaying the exact channel reaches -54 dB in
      % both).  Each threshold sits a few dB above what correct code
      % measures.
      if doppler == 0
        limit = -40;
      else
        limit = -30;
      end
      for start = [2000, late]
        replayed = replay(x, testCase.Fs, 1:M, channel, start);
        for m = 1:M
          [nmse, A, B] = compare(testCase, truth(:, m), replayed(:, m), probe);
          testCase.verifyLessThan(nmse, limit, sprintf(['Element %d, start ', ...
            '%d: replayed output is %.1f dB from the truth.'], m, start, nmse));
          if m == 1 && start == 2000
            testCase.log(1, matlab.unittest.diagnostics.FigureDiagnostic( ...
              plot_round_trip(testCase, A, B, doppler, nmse), 'Formats', 'png'));
          end
        end
      end

      % Control: with the phase convention reversed, the comparison has to
      % fail, or it would not notice a sign error in est.
      if doppler ~= 0
        channel.phi_hat = -phi_hat;
        flipped = replay(x, testCase.Fs, 1, channel, 2000);
        % Aligned as finely as compare now does, the negated phase still
        % leaves the replay far from the truth.
        nmse = compare(testCase, truth(:, 1), flipped, probe);
        testCase.verifyGreaterThan(nmse, -5, ...
          sprintf('With phi_hat negated, the replay is still %.1f dB from the truth.', nmse));
      end
    end

    function testRejectsInvalidInputs(testCase)
      [data, params] = simulate(testCase, 1, 0);

      bad = params;
      bad.optim = 4;
      testCase.verifyError(@() est(data, bad), ?MException);

      bad = params;
      bad.optim = 0;
      testCase.verifyError(@() est(data, bad), ?MException);

      bad = rmfield(params, 'mu');
      bad.optim = 1;
      testCase.verifyError(@() est(data, bad), ?MException);

      bad = rmfield(params, 'K_1');
      testCase.verifyError(@() est(data, bad), ?MException);

      bad = params;
      bad.M = 2;
      testCase.verifyError(@() est(data, bad), ?MException);

      bad = rmfield(data, 'fc');
      testCase.verifyError(@() est(bad, params), ?MException);
    end

    function testPlots(testCase)
      M = 2;
      [data, params] = simulate(testCase, M, 0);
      params.optim = 2;
      [h_hat, phi_hat] = est(data, params);

      channel_params.fs_delay = testCase.R * params.nsd;
      channel_params.fs_time = testCase.R * params.nst;
      channel_params.fc = testCase.fc;

      meta = params;
      meta.description = 'testPlots';
      meta.Fs = testCase.Fs;
      meta.fc = testCase.fc;
      meta.element_spacing = 0.075;
      meta.vertical = true;
      meta.limit = -30;

      figs = plots(h_hat, phi_hat, channel_params, meta, true);
      testCase.addTeardown(@() close(findobj('Type', 'figure')));

      for name = {'main', 'ir', 'doppler', 'mip', 'pll'}
        testCase.verifyTrue(isfield(figs, name{1}) && isgraphics(figs.(name{1})), ...
          sprintf('plots did not return a figure named %s.', name{1}));
      end
    end

  end
end


%% Helpers

function [data, params] = simulate(testCase, M, a, symbols)
% Pass N repetitions of a BPSK sequence through the known channel, with
% Doppler factor a (a > 0: receding), and add noise 40 dB below the signal.
% With SYMBOLS given, pass those once instead, without noise.
N = 6;
noisy = nargin < 4;
if noisy
  d = randi([0, 1], 2^12-1, 1) * 2 - 1;
  symbols = repmat(d, N, 1);
end

s = passband(testCase, symbols);
s = [zeros(testCase.Fs/4, 1); s; zeros(testCase.Fs/4, 1)];
r = zeros(length(s), M);
for m = 1:M
  for p = 1:numel(testCase.delays)
    tau = testCase.delays(p) + (m-1) * testCase.element_delay;
    r(:, m) = r(:, m) + testCase.gains(p) * circshift(s, round(tau * testCase.Fs));
  end
end
[up, down] = rat(1 + a);
r = resample(r, up, down, 'Dimension', 1);
r = r ./ sqrt(mean(pwr(r)));
if noisy
  r = r + 0.01 * randn(size(r));
end

data.Fs = testCase.Fs;
data.R = testCase.R;
data.fc = testCase.fc;
data.r = r;
if noisy
  data.d = d;
  params = default_params(M, N);
end
end


function s = passband(testCase, symbols)
u = resample(symbols, testCase.Fs / testCase.R, 1);
s = real(u .* exp(2j*pi*testCase.fc*(0:length(u)-1).' / testCase.Fs));
end


function s = symbol_train(testCase, symbols, nsd)
% The symbols as est sees them, a zero-stuffed train at fs_delay = nsd*R,
% carried to Fs and up to fc, so that replay applies h_hat to them once.
x = zeros(numel(symbols) * nsd, 1);
x(1:nsd:end) = symbols;
u = resample(x, testCase.Fs / (testCase.R * nsd), 1);
s = real(u .* exp(2j*pi*testCase.fc*(0:length(u)-1).' / testCase.Fs));
end


function params = default_params(M, N)
params.nsd = 2;
params.nst = 2;
params.nslr = 2;
params.M = M;
params.N = N - 1;
params.K_1 = 10;
params.K_2 = 50;
params.delta = 0;
params.sync_length = 64;
params.delay_tracking = true;
params.optim = 2;
params.mu = 0.3 / (params.K_1 + params.K_2) / params.nsd;
params.lambda = 0.995;
params.regularization = 1e-3;
params.Kf_1 = 1e-3;
params.Kf_2 = params.Kf_1 / 10;
end


function h = expected_response(testCase, params, offset)
% The baseband impulse response the simulated channel produces, sampled on
% est's tap grid: the transmit pulse at each path delay, times the path gain
% and the carrier phase exp(-j*2*pi*fc*tau).  Tap K_1*nsd + 1 + OFFSET is
% delay zero.
fs_delay = testCase.R * params.nsd;
[q, tq] = transmit_pulse(testCase);
k = (1:(params.K_1 + params.K_2) * params.nsd).';
t = (k - params.K_1*params.nsd - 1 - offset) / fs_delay;
h = zeros(numel(k), params.M);
for m = 1:params.M
  for p = 1:numel(testCase.delays)
    tau = testCase.delays(p) + (m-1) * testCase.element_delay;
    h(:, m) = h(:, m) + testCase.gains(p) * exp(-2j*pi*testCase.fc*tau) ...
      * interp1(tq, q, t - tau, 'spline', 0);
  end
end
end


function [nmse, offset] = fit_response(testCase, h, params, conjugate)
% Best fit of the expected response to H over one tap offset shared by all
% elements and one complex gain per element, as normalized squared error in
% dB.  Each element has its own PLL, which settles to its own constant
% phase, so the gain cannot be shared.
nmse = inf;
offset = NaN;
for o = -2:0.02:2
  ref = expected_response(testCase, params, o);
  if conjugate
    ref = conj(ref);
  end
  g = sum(conj(ref) .* h, 1) ./ sum(abs(ref).^2, 1);
  e = sum(abs(h - g .* ref).^2, 'all') / sum(abs(h).^2, 'all');
  if e < nmse
    nmse = e;
    offset = o;
  end
end
nmse = 10*log10(nmse);
end


function [q, t] = transmit_pulse(testCase)
% The pulse simulate transmits for one symbol, and its time axis with zero
% at the symbol.
n = 20;
impulse = zeros(2*n + 1, 1);
impulse(n + 1) = 1;
up = testCase.Fs / testCase.R;
q = resample(impulse, up, 1);
t = ((0:numel(q)-1).' - n*up) / testCase.Fs;
end


function tap = peak_tap(power)
% Power-weighted position of the strongest peak, to a fraction of a tap.
[~, k] = max(power);
w = max(k-2, 1):min(k+2, numel(power));
tap = sum(w(:) .* power(w(:))) / sum(power(w));
end


function [nmse, A, B] = compare(testCase, truth, replayed, probe)
% Align the replayed output to the truth on the probe's correlation peak,
% then measure the error at baseband, after the best complex gain, in dB.
% A and B are the aligned truth and replay at 4R, B scaled by that gain.
Fs = testCase.Fs;
template = resample(probe(1:64), Fs / testCase.R, 1);
v_truth = truth .* exp(-2j*pi*testCase.fc*(0:numel(truth)-1).' / Fs);
v_replay = replayed .* exp(-2j*pi*testCase.fc*(0:numel(replayed)-1).' / Fs);
[x1, l1] = xcorr(v_truth, template);
[x2, l2] = xcorr(v_replay, template);
[~, k1] = max(abs(x1));
[~, k2] = max(abs(x2));
start_truth = l1(k1) + 1;
start_replay = l2(k2) + 1;

n = 3 * 2^14;
A = lowpass_bb(v_truth(start_truth + (0:n-1)), testCase);
B = lowpass_bb(v_replay(start_replay + (0:n-1)), testCase);

% The two starts can differ by a fraction of a sample.  Under Doppler a
% search in eighths of a sample leaves up to 1.3 us of timing error, which
% alone reads about -35 dB, so refine the best eighth continuously.
k = (0:numel(B)-1).';
err = @(shift) fit_error(A, interp1(k, B, k + shift, 'spline', 0));
coarse = -1:0.125:1;
[~, i] = min(arrayfun(err, coarse));
shift = fminbnd(err, coarse(max(i-1, 1)), coarse(min(i+1, numel(coarse))), ...
  optimset('TolX', 1e-6));
[e, B] = err(shift);
nmse = 10*log10(e);
end


function [e, gB] = fit_error(A, B)
g = (B' * A) / (B' * B);
gB = g * B;
e = sum(abs(A(100:end-100) - gB(100:end-100)).^2) / sum(abs(A(100:end-100)).^2);
end


function fig = plot_round_trip(testCase, A, B, doppler, nmse)
% The truth and the replay at baseband over a few milliseconds, and the
% error along the whole comparison.
fs = 4 * testCase.R;
t = (0:numel(A)-1).' / fs * 1e3;
fig = figure;
subplot(2, 1, 1);
w = t >= 100 & t < 103;
plot(t(w), real(A(w)), 'LineWidth', 1.2), hold on
plot(t(w), real(B(w)), '--', 'LineWidth', 1.2)
xlabel('Time [ms]'), ylabel('Re\{v\}')
legend('Truth', 'Replay of the estimate', 'Location', 'best')
title(sprintf('Round trip, Doppler %+.0e: NMSE %.1f dB', doppler, nmse))
grid on
subplot(2, 1, 2);
e = abs(A - B).^2 / mean(abs(A).^2);
plot(t, 10*log10(movmean(e, 256)))
xlabel('Time [ms]'), ylabel('Error [dB]')
title('Squared error relative to the mean signal power, moving average')
grid on
end


function v = lowpass_bb(v, testCase)
% Resample to 4R.  The anti-aliasing filter keeps the signal band and drops
% the image at -2*fc that downshifting the passband leaves.
[up, down] = rat(4 * testCase.R / testCase.Fs);
v = resample(v, up, down);
end


function p = pwr(x)
p = mean(abs(x).^2, 1);
end
