# Fast S-box algebraic degree computation

Code and experimental data for **A Series of Faster Approaches for Computing the Algebraic Degree Properties of S-boxes**, by Renjie Zhou, Zhen Li and Yan Tong.

The implementations compute the minimum algebraic degree and update the maximum algebraic degree or the complete component-degree spectrum after an output transposition. The artifact includes single-thread experiments, complete 8-bit local searches and CPU parallel scaling.

## Start here

On Linux x86-64 with GCC 13 or later, GNU Make and Python 3.10 or later:

```sh
git clone https://github.com/MarstonRaines/fast-sbox-degree.git
cd fast-sbox-degree
make -j4 all
make check
python3 scripts/analyze.py --check
```

The build downloads 21 unchanged PEIGEN headers from a fixed upstream commit and checks every file against a saved SHA-256 digest. `make check` runs independent small-input checks, update/rollback checks, exhaustive 2-by-2 mappings, component enumeration checks and forced parallel execution. `make check-full` additionally checks every stored component-ANF coefficient in 36 states of two 16-bit inputs; it uses up to 64 threads.

The 16-bit spectrum state needs substantial memory: about 512 MiB with packed coefficients, versus about 16 GiB for the scalar reference representation, before auxiliary buffers. Full reference timing sweeps take hours. Reading the existing results and recomputing their tables needs no C++ build.

## Contents

| Directory | Contents |
|---|---|
| `src/serial/` | Scalar reference algorithms, packed proposed algorithms and the PEIGEN adapter |
| `src/parallel/` | OpenMP implementation of the same numerical tasks |
| `data/inputs/` | 147 explicit LUTs and fixed swap/check sequences |
| `data/measurements/` | 17,170 timing records and 1,460 separate operation-count records |
| `data/outputs/` | Deduplicated complete final LUTs, degree histograms and accepted swaps |
| `protocol/` | Exact ordered tasks, repetitions, thread caps and frozen workload sizes |
| `results/serial/` | Reference, optimized, practical-instance and search tables |
| `results/parallel/` | Thread scaling, equal-thread algorithm comparisons and serial anchors |
| `scripts/` | Dependency download, workload replay, table reconstruction and plotting |
| `figures/` | Eight publication figures and one numerical table in PDF, SVG and 300 dpi PNG |
| `docs/` | Protocol, results, measured-build provenance and release validation |

See [the experimental protocol](docs/protocol.md) for the distinction between single-evaluation, fixed-update and full-search workloads, and [the results guide](docs/results.md) for the interpretation of speedups.

## Reproduce measurements

Start with a small subset:

```sh
python3 scripts/run.py --suite serial --group B --case random_n08_00 --metric spectrum --limit 2 --output runs/smoke
python3 scripts/run.py --suite parallel --group B --case random_n16_00 --metric spectrum --threads 4 --limit 2 --output runs/parallel-smoke
```

Use `--list` to inspect selected tasks; omit filters to replay a full suite. Optional `--cpu-start` and `--numa-node` use `numactl` for binding. Select consecutive physical cores from your own machine's topology. Every run uses a new output directory and checks its full final output against the saved expected output.

Regenerate the published result tables:

```sh
python3 scripts/analyze.py
```

These tables are derived from the archived measurements in `data/measurements/`. New timings are written only to the chosen `runs/` directory and do not replace the archived measurements.

## Figures

```sh
python3 -m pip install -r requirements.txt
python3 scripts/plot_results.py
```

The [figure guide](docs/figures.md) contains the current previews, source-table links, statistical definitions and complete captions. Each metric has a time/speedup pair with direct numerical labels. Parallel figures compare the methods at equal thread limits; within-method scaling is tabulated separately. All plots read the saved result tables directly.

![Minimum-degree time and speedup](figures/01_minimum_degree.png)

![Complete search time](figures/04_search_time.png)

## Version and license

Release **v1.0.1** updates figure presentation: separate time/speedup pairs, direct numerical labels, search-time bars, equal-thread parallel comparisons and a practical-instance table. The implementation, protocols and archived measurements are unchanged from v1.0.0.

This is a publication packaging revision: the retained numerical kernels match the measured source versions, while validation and reporting have been repackaged. [Provenance](docs/provenance.json) records the measured source and binary digests; [release validation](docs/release-validation.json) records checks of the packaged build. The archived timings were not regenerated using the packaged build.

The project is licensed under **GPL-3.0-only**; see [LICENSE](LICENSE). PEIGEN is an external dependency with its own retained [GPL-3.0 license](third_party/PEIGEN-LICENSE), authorship and header notices. Its exact revision is recorded in [third_party/peigen.json](third_party/peigen.json). Please also cite Bao, Guo, Ling and Sasaki's [PEIGEN paper](https://eprint.iacr.org/2019/209) when using that baseline.
