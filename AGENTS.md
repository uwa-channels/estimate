# AGENTS.md

Notes for coding agents working in this repository.  Read `README.md` for what
the scripts do and for the output format.

## Layout

- `src/est.m` estimates the channel.  Its help block documents every field of
  `data` and `params`.
- `src/plots.m` plots an estimated or a saved channel.
- `examples/example_estimate.m` is the end-to-end example: it makes a
  synthetic channel, estimates it, saves it in the library format, and plots
  it.
- `install.m` adds `src/` and `examples/` to the search path and saves it.

## Running it

```bash
matlab -batch "addpath('src', 'examples'); example_estimate"
octave --eval "pkg load signal statistics; addpath('src', 'examples'); example_estimate"
```

There is no test suite.  Run `example_estimate.m` after a change, and check that the
estimated impulse response and PLL phase in the plots look right.  It needs the
Signal Processing Toolbox in MATLAB (R2021a or later), or the `signal` and
`statistics` packages in Octave 9.

## File format

`example_estimate.m` saves channels in the format that `uwa-channels/matlab` and
`uwa-channels/python` read.  The specification is in `CONTRIBUTING.md` of
those repositories.  When `params.delay_tracking` is true, save the phase
output of `est` as `phi_hat`; otherwise save it as `theta_hat`.  Do not change
the saved fields or their axis order unless both toolboxes change with them.

## Conventions

- Each function starts with a help block: an H1 line, the call forms, Inputs,
  Outputs, Example, required files and toolboxes, See also, Author, Email,
  License, and a dated Revision history.
- Two-space indentation.
- Commit messages are one imperative sentence that says what changed and why.

## Citation

`CITATION.cff` cites the arXiv paper (arXiv:2609.03207) and the Zenodo data
(doi:10.5281/zenodo.21287414).
