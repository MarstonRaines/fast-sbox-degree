# Fast S-box algebraic degree computation

Code and experimental data for **A Series of Faster Approaches for Computing the Algebraic Degree Properties of S-boxes**, by Renjie Zhou, Zhen Li and Yan Tong.

The implementations compute the minimum algebraic degree and update the maximum algebraic degree or the full degree spectrum after an output transposition. All three metrics are evaluated at every width from 3 to 19 bits. A second experiment improves three published 8-bit S-boxes from minimum degree 6 to 7 while preserving standard vectorial nonlinearity exactly at 104, and compares total target times with PEIGEN.

## Start here

On Linux x86-64 with GCC 13 or later, GNU Make and Python 3.10 or later:

```sh
git clone https://github.com/MarstonRaines/fast-sbox-degree.git
cd fast-sbox-degree
make -j4 all
make check
python3 scripts/analyze.py --check
```

The build downloads 21 unchanged PEIGEN headers from a fixed upstream commit and checks every file against a saved SHA-256 digest. `make check` runs independent small-input checks, update/rollback checks, exhaustive 2-by-2 mappings and component-enumeration checks. `make check-full` additionally checks every stored component-ANF coefficient in 12 states of two 16-bit inputs.

The packed spectrum state uses about 512 MiB at 16 bits and 32 GiB at 19 bits, before auxiliary buffers. Minimum/maximum computations use much smaller states. Reading the archived results and recomputing their tables needs no C++ build.

## Contents

| Directory | Contents |
|---|---|
| `src/serial/` | Preserved implementation for the reference timings and counts |
| `src/optimized/`, `src/extended/`, `src/minmax/` | Frozen measured implementations for 3–16 bits, 17–19-bit spectrum, and 17–19-bit minimum/maximum, respectively |
| `src/postprocess.cpp` | Common strict-NL literature postprocessing driver |
| `data/inputs/` | 177 explicit LUTs: 170 random, six practical, one legacy identity control |
| `data/literature/` | Three published starting LUTs, improved LUTs, source attribution and accepted trajectories |
| `data/measurements/` | 13,000 reference/fixed-workload/legacy-search timings, 30 literature timings, and 1,460 separate counts |
| `data/outputs/` | Deduplicated complete final LUTs, degree histograms and accepted swaps |
| `protocol/` | Exact ordered tasks, repetitions and frozen workload sizes |
| `results/serial/` | Reference, optimized, practical-instance and search tables |
| `scripts/` | Dependency download, workload replay, table reconstruction and plotting |
| `figures/` | Four publication figures and one numerical table in PDF, SVG and 300 dpi PNG |
| `docs/` | Protocol, results, measured-build provenance and release validation |

See [the experimental protocol](docs/protocol.md) for the distinction between single-evaluation, fixed-update and full-search workloads, and [the results guide](docs/results.md) for the interpretation of speedups.

## Reproduce measurements

Start with a small subset:

```sh
python3 scripts/run.py --suite serial --group B --case random_n08_00 --metric spectrum --limit 2 --output runs/smoke
python3 scripts/literature.py --check
python3 scripts/literature.py --replay --output runs/literature
```

Use `--list` to inspect selected tasks; omit filters to replay the full serial suite. Optional `--cpu-start` and `--numa-node` use `numactl` for binding. Every run uses a new output directory and checks its full final output against the saved expected output.

The literature checker independently recomputes full-component NL, degree spectra, differential uniformity, and the degree-7 coefficient ranks along every accepted path. `--replay` additionally reruns both C++ evaluators and checks their complete candidate digests and final LUTs. The shared nonlinearity checks are included in total search times.

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

The [figure guide](docs/figures.md) contains the current previews, source-table links, statistical definitions and complete captions. Each metric has a time/speedup pair with direct numerical labels. All plots read the saved result tables directly.

![Minimum-degree time and speedup](figures/01_minimum_degree.png)

![Complete search time](figures/04_search_time.png)

## Version and license

Release **v1.1.0** contains the optimized 3–19-bit comparison for all three metrics, six practical components, and three literature postprocessing cases at exact NL 104. The earlier unconstrained 14-start searches remain as additional archived experiments; Figure 4 now presents the literature target times. Older releases and their measurements remain available in Git history and release assets.

The separate source directories preserve the measured numerical versions. [Provenance](docs/provenance.json) identifies their roles and source/binary digests; [release validation](docs/release-validation.json) records checks of the packaged build. Replays check numerical identity and produce new timings without overwriting the archive.

The project is licensed under **GPL-3.0-only**; see [LICENSE](LICENSE). PEIGEN is an external dependency with its own retained [GPL-3.0 license](third_party/PEIGEN-LICENSE), authorship and header notices. Its exact revision is recorded in [third_party/peigen.json](third_party/peigen.json). Please also cite Bao, Guo, Ling and Sasaki's [PEIGEN paper](https://doi.org/10.13154/tosc.v2019.i1.330-394) when using that baseline.
