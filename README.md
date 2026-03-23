# Underwater Acoustic Channel Toolbox — Channel Estimation Scripts

[![Generic badge](https://img.shields.io/badge/MATLAB-R2021a-BLUE.svg)](https://shields.io/) with [Signal Processing Toolbox™](https://www.mathworks.com/products/signal.html) or [![Generic badge](https://img.shields.io/badge/Octave-9.0-BLUE.svg)](https://shields.io) with the [signal](https://gnu-octave.github.io/packages/signal/) and [statistics](https://gnu-octave.github.io/packages/statistics/) packages.

This MATLAB®/Octave toolbox estimates underwater acoustic channels from single-carrier recordings and visualizes the results.

Please report bugs and suggest enhancements by [creating a new issue](https://github.com/uwa-channels/estimate/issues). We welcome your feedback.

## Quick start

Run `main.m` in MATLAB or Octave.  The script generates a BPSK signal, passes it through a synthetic multipath channel with Doppler, estimates the time-varying impulse response, and plots the results.

## Functions

| Function | Description |
|----------|-------------|
| `est.m` | Adaptive channel estimation with LMS, RLS, or SFTF, and optional delay tracking via a second-order PLL. |
| `plots.m` | Visualization of the estimated channel: time-varying IR, Doppler-delay scattering function, delay-angle (broadband beamforming), multipath intensity profile, PLL phase, and estimation parameters. |

## Output format

The estimator produces two structs that are saved alongside the channel data:

**`params`** — physical parameters of the channel.

| Field | Description |
|-------|-------------|
| `fs_delay` | Sampling rate in the delay domain [Hz] |
| `fs_time` | Sampling rate along the time axis [Hz] |
| `fc` | Center frequency [Hz] |

**`meta`** — estimation metadata.

| Field | Description |
|-------|-------------|
| `description` | Free-text description of the experiment |
| `nsd` | Samples per symbol in the delay domain |
| `nst` | Samples per symbol in the time domain |
| `Fs`  | Sample rate of the received signal |
| `K_1`, `K_2` | Anti-causal and causal filter lengths [symbols] |
| `fc` | Center frequency [Hz] |
| `element_spacing` | Array element spacing [m] |
| `vertical` | `true` if vertical array |
| `delay_tracking` | `true` if delay tracking is enabled |
| `limit` | Lower dB limit for plotting |
| `optim` | Optimizer (1: LMS, 2: RLS, 3: SFTF) |
| `mu` | LMS step size (when `optim == 1`) |
| `lambda` | Forgetting factor (when `optim == 2` or `3`) |
| `Kf_1`, `Kf_2` | PLL loop filter coefficients |
| `nslr` | Delay tracking rate (when `delay_tracking == true`) |
| `codename` | Short identifier for the channel |

The `.mat` file contains `h_hat`, either `theta_hat` or `phi_hat` (depending on the tracking mode), `params`, `meta`, and `version`.

## Plotting

`plots.m` returns a struct of figure handles.  Saving is left to the caller.

```matlab
% After estimation:
figs = plots(h_hat, theta_hat, params, meta);

% With standalone subplot figures:
figs = plots(h_hat, theta_hat, params, meta, true);

% Export:
exportgraphics(figs.main, "channel.pdf", "Resolution", 300);
exportgraphics(figs.ir, "channel_ir.pdf", "Resolution", 300);
```

To plot a previously saved channel:

```matlab
channel = load('yellow_2.mat');
if isfield(channel, 'phi_hat')
  phase = channel.phi_hat;
else
  phase = channel.theta_hat;
end
figs = plots(channel.h_hat, phase, channel.params, channel.meta);
```

The returned struct `figs` may contain the following fields:

| Field | Description |
|-------|-------------|
| `main` | 2×3 summary figure |
| `ir` | Time-varying impulse response |
| `doppler` | Doppler-delay scattering function |
| `delay_angle` | Delay-angle plot (vertical arrays with M ≥ 13) |
| `mip` | Multipath intensity profile |
| `pll` | PLL phase trajectory |


## License
The license is available in the [LICENSE](LICENSE) file within this repository.

© 2025-2026, Underwater Acoustic Channels Group.
