# Underwater Acoustic Channel Toolbox — Channel Estimation Scripts

[![Generic badge](https://img.shields.io/badge/MATLAB-R2021a-BLUE.svg)](https://shields.io/) with [Signal Processing Toolbox™](https://www.mathworks.com/products/signal.html) or [![Generic badge](https://img.shields.io/badge/Octave-9.0-BLUE.svg)](https://shields.io) with the [signal](https://gnu-octave.github.io/packages/signal/) and [statistics](https://gnu-octave.github.io/packages/statistics/) packages.

This MATLAB®/Octave toolbox estimates an underwater acoustic channel using single-carrier signals.

Please report bugs and suggest enhancements by [creating a new issue](https://github.com/uwa-channels/estimate/issues). We welcome your feedback.

## Using the channel estimation script

Simply run `estimate.m` in MATLAB or Octave. The script generates a single-carrier BPSK-modulated signal, passes it through a multipath channel, and then, on the receive side, estimates the time-varying channel impulse responses. Finally, it produces an image showing the magnitude of the estimated time-varying channel impulse responses.

The channel estimation and delay tracking mechanism is described [here](https://ieeexplore.ieee.org/abstract/document/10942614).
