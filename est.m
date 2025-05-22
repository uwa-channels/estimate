function [h_hat_store, theta_hat_store] = est(data, params)
%EST Estimate channel.
%   Detailed explanation TBD.
%
% Author: Zhengnan Li
% Email : uacr@ofdm.link
%
% License: MIT
%

%% Check inputs
validate_inputs(data, params);

%% Unpack data and parameters
Fs = data.Fs; % Original data sampling rate [samples per second]
R = data.R; % Symbol rate [Baud]
fc = data.fc; % Center frequency [Hz]
d = data.d; % Data symbols, num_symbols x 1
r = data.r; % Received data symbols, num_channels x num_samples

nsd = params.nsd; % Number of samples per symbol, in delay domain
nst = params.nst; % Number of samples per symbol, in time domain
M = params.M; % Number of receiving channels
N = params.N; % Number of data repeats
delay_tracking = params.delay_tracking;

K_1 = params.K_1; % Number of symbols to the left (anti-causal)
K_2 = params.K_2; % Number of symbols to the right (causal)
K = K_1 + K_2; % Total number of symbols
L = K * nsd; % Filter length

delta = params.delta; % Synchronization fine adjustment
Ts = 1 / (R * nsd); % Sample interval, delay [second]
Tr = Ts / params.nslr;

sync_length = params.sync_length;

if params.delay_tracking
    fs = 1 / Tr;
else
    fs = 1 / Ts;
end

col = (1:M).'; % Dummy variable for finding indices

optim = params.optim;
if optim == 1 % LMS
    mu = params.mu;
elseif optim == 2 % RLS
    lambda = params.lambda;
    P = (1 / params.regularization) * eye(nsd*(K_1 + K_2));
    P = repmat(P, 1, 1, M);
    errors = zeros(1, 1, M);
    inputs = zeros(K*nsd, 1, M);
elseif optim == 3 % SFTF
    lambda = params.lambda;
    inv_delta = params.regularization;
    gaminv = params.conversion * ones(1, M);
    w = zeros(L, M);
    uN = zeros(L, M);
    wf = zeros(L, M); % forward prediction vector wf
    wb = zeros(L, M); % backward prediction vector wb
    g = zeros(L, M); % gain vector
    xifinv = lambda / inv_delta * ones(1, M); % forward prediction energy
    xib = inv_delta / lambda^(L + 1) * ones(1, M); % backward prediction energy
    k = [1.5, 2.5, 1]; % gains for convex combinations
end

Kf_1 = params.Kf_1; % PLL coefficients 1
Kf_2 = params.Kf_2; % PLL coefficients 2

%% Downshift, resample, find correlation peaks, and align the signal
% Downshift
v_bb = r .* exp(-1j*2*pi*fc*(0:size(r, 1) - 1).'./Fs);

% Find correlation peaks
[p, q] = rat(Fs/R);
[xcor, lags] = xcorr(v_bb(:, 1), resample(d(1:sync_length), p, q));
xcor = abs(xcor) / max(abs(xcor));
xcor(lags < 0) = [];
[~, locs] = findpeaks(xcor, 'MinPeakHeight', 0.3, ...
    'MinPeakDistance', length(d)*(Fs / R)-100);
% lags(lags<0) = [];
% figure, hold on, plot(lags, xcor)
% plot(lags(locs), xcor(locs), 'x');

% Align the signal
v_bb = v_bb(locs(1)+delta:end, :);

% Resample
[up, down] = rat(fs/Fs);
v = resample(v_bb, up, down, 'Dimension', 1);
size_v = size(v);

% Normalize the signal
powers = sum(abs(v).^2, 1) ./ size(v, 1);
v = v ./ sqrt(sum(powers));

%% Channel estimation loop
d_all = zeros(N*length(d)*nsd, 1);
% Repeat the data symbols
d_all(1:nsd:end) = repmat(d, N, 1);

% Pre-allocation
h_hat_store = zeros(K*nsd, M, floor((length(d) * N - K)*nst));
h_hat = zeros(K*nsd, M);
theta_hat = zeros(M, length(d_all)-K_1*nsd);
sum_psi = zeros(M, 1);

% Loop over the data symbols
count = 1;
error = zeros(length(d_all)-K*nsd, M);
for n = K_2 * nsd:length(d_all) - K_1 * nsd
    if mod(n/nsd, 10000) == 0
        fprintf("n = %d\n", n/nsd);
    end

    % Save the data every 1/nst samples
    if mod(n/nsd, 1/nst) == 0
        h_hat_store(:, :, count) = h_hat;
        count = count + 1;
    end

    window = -K_2 * nsd + (K * nsd:-1:1) + n;
    d_vec = d_all(window);

    if ~delay_tracking % Without delay tracking
        v_hat = h_hat' * d_vec;
        err = v(n, :).' .* exp(-1j*theta_hat(:, n)) - v_hat;
    else % With delay tracking
        v_hat = h_hat' * d_vec;
        tn = n * Ts - theta_hat(:, n) ./ (2 * pi * fc);
        tln = floor(tn/Tr) * Tr;
        trn = ceil(tn/Tr) * Tr;
        alpha = (tn - tln) / Tr;
        y = (1 - alpha) .* v(sub2ind(size_v, floor(tln/Tr), col)) + ...
            alpha .* v(sub2ind(size_v, ceil(trn/Tr), col));
        err = y .* exp(-1j*theta_hat(:, n)) - v_hat;
    end
    error(n-K_2*nsd+1, :) = err;

    % Adaptation
    if optim == 1
        h_hat = h_hat + update_lms(d_vec, err, mu);
    elseif optim == 2
        errors(1, 1, :) = err;
        inputs(:, 1, :) = d_vec .* ones(1, M);
        [W, P] = update_rls(inputs, errors, lambda, P);
        h_hat = h_hat + W;
    elseif optim == 3
        inp = d_vec(1) .* ones(1, M); % Input
        fa = inp - sum(uN.*wf); % a-priori forward error
        f = fa ./ gaminv; % a-posteriori forward error
        % Order and time update of gain vector
        g1 = [zeros(1, M); g] + conj(fa) .* xifinv ./ lambda .* [ones(1, M); -wf];
        v1 = g1(L+1, :);
        v2 = g1(1, :);
        gaminva = gaminv + v2 .* conj(fa);
        xifinv = xifinv ./ lambda - abs(v2).^2 ./ gaminva;
        wf = wf + f .* g;
        ba1 = lambda * xib .* v1;
        ba2 = uN(L, :) - sum([inp; uN(1:L-1, :)].*wb);
        ba3_1 = ba2 * k(1) + ba1 * (1 - k(1));
        ba3_2 = ba2 * k(2) + ba1 * (1 - k(2));
        ba3_3 = ba2 * k(3) + ba1 * (1 - k(3));
        gaminv1 = gaminva - v1 .* conj(ba3_3);
        b_1 = ba3_1 ./ gaminv1;
        b_2 = ba3_2 ./ gaminv1;
        xib = lambda * xib + b_2 .* conj(ba3_2);
        g = g1(1:L, :) + v1 .* wb;
        wb = wb + b_1 .* g;
        gaminv = 1 ./ (lambda^(L) * xifinv .* xib);
        uN = [inp; uN(1:L-1, :)];
        ep = err.' ./ gaminv; % a-posteriori estimation error
        w = w + conj(ep) .* g; % time-update of weight vector
        h_hat = w;
    end

    % PLL updates
    psi = -imag(v_hat.*conj(err));
    sum_psi = sum_psi + psi;
    theta_hat(:, n+1) = theta_hat(:, n) + Kf_1 .* psi + Kf_2 .* sum_psi;
end

theta_hat_store = theta_hat;
if optim == 1
    h_hat_store = h_hat_store(:, :, ceil(20*K*nst):end); % LMS converges after 20 * K
    theta_hat_store = theta_hat(:, ceil(20*K*nsd):end);
    h_hat_store = conj(h_hat_store);
elseif optim == 2
    h_hat_store = h_hat_store(:, :, ceil(2*K*nst):end); % RLS converges sooner
    theta_hat_store = theta_hat(:, ceil(2*K*nsd):end);
    h_hat_store = conj(h_hat_store);
elseif optim == 3
    h_hat_store = h_hat_store(:, :, ceil(2*K*nst):end); % SFTF converges sooner
    theta_hat_store = theta_hat(:, ceil(2*K*nsd):end);
    h_hat_store = conj(h_hat_store);
end

% Truncate the delay and phase
h_hat_store = h_hat_store(:, :, 1:end-1); % Ignore last sample due to phase
theta_hat_store = theta_hat_store(:, 1:floor(size(h_hat_store, 3)*nsd/nst));

end

%% Update functions
function W = update_lms(input, err, mu)
W = mu * input .* err';
end

function [W, P] = update_rls(inputs, errors, lambda, P)
common = (1 / lambda) * pagemtimes(P, inputs);
kal = common ./ (1 + pagemtimes(pagectranspose(inputs), common));
P = (1 / lambda) * P - (1 / lambda) * ...
    pagemtimes(kal, pagemtimes(pagectranspose(inputs), P));
W = squeeze(pagemtimes(kal, pagectranspose(errors)));
end

%% Validate inputs
function validate_inputs(data, params)

% Check data
assert(all(isfield(data, {'Fs', 'R', 'fc', 'd', 'r'})), ...
    'The struct "data" does not have all required fields.');

% Check params
assert(all(isfield(params, {'nsd', 'nst', 'M', 'N', ...
    'K_1', 'K_2', 'Kf_1', 'Kf_2', 'delta', 'delay_tracking'})), ...
    'The struct "params" does not have all required fields.');

if params.optim == 1
    assert(all(isfield(params, {'mu'})), ...
        'Optimizer 1 does not have all required parameters.');
    convergence_period = 20 * (params.K_1 + params.K_2) * params.nsd;
elseif params.optim == 2
    assert(all(isfield(params, {'lambda', 'regularization'})), ...
        'Optimizer 2 does not have all required parameters.');
    convergence_period = 2 * (params.K_1 + params.K_2) * params.nsd;
elseif params.optim == 3
    assert(all(isfield(params, {'lambda', 'regularization', 'conversion'})), ...
        'Optimizer 3 does not have all required parameters.');
    convergence_period = 2 * (params.K_1 + params.K_2) * params.nsd;
elseif params.optim == 4
    assert(all(isfield(params, {'lambda_ch'})), ...
        'Optimizer 4 does not have all required parameters.');
    convergence_period = 20 * (params.K_1 + params.K_2) * params.nsd;
else
    error('Invalid optimizer option.');
end


% Validate data
validateattributes(data.Fs, 'numeric', {'scalar', 'positive', 'finite'}, ...
    '', 'data.Fs (signal sampling rate [samples per second])');
validateattributes(data.R, 'numeric', {'scalar', 'positive', 'finite'}, ...
    '', 'data.R (signal symbol rate [Baud])');
validateattributes(data.fc, 'numeric', {'scalar', 'positive', 'finite'}, ...
    '', 'data.fc (signal center frequency [Hz])');
validateattributes(data.d, 'numeric', {'ncols', 1, 'finite'}, ...
    '', 'data.d (data vector, n_data x 1)');
validateattributes(data.r, 'numeric', {'ncols', params.M}, ...
    '', 'data.r (signal samples, n_samples x number of receivers)');

% Validate params
validateattributes(params.N, 'numeric', {'scalar', 'positive', 'integer', ...
    '<=', floor(size(data.r, 1)/length(data.d)/(data.Fs / data.R)), ...
    '>=', ceil((convergence_period + params.K_1 * params.nsd)/length(data.d))}, ...
    '', 'params.N (number of repeatations of the data symbols)');
validateattributes(params.nsd, 'numeric', {'scalar', '>=', 1, '<=', data.Fs / data.R}, ...
    '', 'params.nsd (number of samples per symbol in delay domain)');
validateattributes(params.nst, 'numeric', {'scalar', '>=', 1 / (length(data.d) * params.N), '<=', params.nsd}, ...
    '', 'params.nst (number of samples per symbol in time domain)');
validateattributes(params.M, 'numeric', {'scalar', 'positive', 'finite'}, ...
    '', 'params.M (number of receivers)');
validateattributes(params.K_1, 'numeric', {'scalar', 'nonnegative', 'finite'}, ...
    '', 'params.K_1 (length of the channel estimator, anti-causal)');
validateattributes(params.K_2, 'numeric', {'scalar', 'nonnegative', 'finite'}, ...
    '', 'params.K_2 (length of the channel estimator, causal');
validateattributes(params.delta, 'numeric', {'scalar', 'finite'}, ...
    '', 'params.delta (fractionasl spacing adjustment)');
if params.optim == 1
    validateattributes(params.mu, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.mu (LMS update factor)');
elseif params.optim == 2
    validateattributes(params.lambda, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.lambda (RLS forgetting factor)');
    validateattributes(params.regularization, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.lambda (RLS regularization factor)');
elseif params.optim == 3
    validateattributes(params.lambda, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.lambda (SFTF forgetting factor)');
    validateattributes(params.regularization, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.lambda (SFTF regularization factor)');
    validateattributes(params.conversion, 'numeric', {'scalar', 'positive', 'finite'}, ...
        '', 'params.lambda (SFTF conversion factor)');
else
    error('Invalid params.optim.')
end
validateattributes(params.Kf_1, 'numeric', {'scalar', 'finite'}, ...
    '', 'params.Kf_1 (PLL coefficient)');
validateattributes(params.Kf_2, 'numeric', {'scalar', 'finite'}, ...
    '', 'params.Kf_2 (PLL coefficient)');
validateattributes(params.delay_tracking, 'logical', {}, '', 'params.delay_tracking');

end

% [EOF]
