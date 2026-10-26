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
- `tests/testEst.m` holds the `matlab.unittest` tests.
- `install.m` adds `src/` and `examples/` to the search path and saves it.

## Running it

```bash
matlab -batch "addpath('src', 'examples'); example_estimate"
octave --eval "pkg load signal statistics; addpath('src', 'examples'); example_estimate"
```

## Tests

```bash
matlab -batch "addpath('src', 'tests', '../matlab/src'); assertSuccess(runtests('tests'))"
```

`testReplayRoundTrip` needs `replay.m` from `uwa-channels/matlab`; the command
above assumes a clone named `matlab` next to this one.  Without it, that test
is skipped.  CI (`.github/workflows/ci.yml`) runs the suite on R2021a and on
the latest release, so keep the code working on R2021a.  The tests do not run
in Octave, which lacks `matlab.unittest`.

The tests compare `h_hat` with the impulse response the simulated channel
must produce.  If a change makes them fail, find out why before you change a
threshold: each one sits a few dB from what correct code measures, and the
header of `testEst.m` lists what each test checks.

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
