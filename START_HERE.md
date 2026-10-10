# Start here

This archive contains the complete DynaComp 0.1.0 source, examples, mathematical tests, retained receipts, GitHub Actions configuration and monograph. Unzip it, open a terminal in `dynamics-compiler`, and run `python -m pip install -e .`. On Windows the Python launcher may be `py`; use `py -m pip install -e .` and substitute `py` for the example commands if necessary.

Run `python -m unittest discover -s tests -v` to exercise the mathematical regressions, then `python examples/run_examples.py` to regenerate the deterministic receipts. Read `README.md` and `docs/MONOGRAPH.md`. The optional performance demonstration needs NumPy and runs with `python examples/benchmark_ring.py`; its timing receipt is machine-dependent.

The repository is intended to be published as a separate project named `dynamics-compiler`. It is independent of PerfectPower. The archive does not include Git credentials or a hidden dependency on the original conversations. All exact mechanisms use the standard library, and the current mathematical limitations are recorded in `docs/DEVELOPMENT.md`.
