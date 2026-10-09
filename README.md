# Fast S-box algebraic degree computation

Code, experimental data and examples for *A Series of Faster Approaches for Computing the Algebraic Degree Properties of S-boxes*, by Renjie Zhou, Zhen Li and Yan Tong.

The implementation computes minimum component degree and complete degree spectra, and updates maximum degree and degree spectra after transpositions. Comparisons use PEIGEN at 3–8 bits and bitwise recomputation at 9–19 bits. The search examples preserve nonlinearity 104 and differential uniformity 8.

## Getting started

```sh
git clone https://github.com/MarstonRaines/fast-sbox-degree.git
cd fast-sbox-degree
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Python 3.11 or later is required for analysis and plotting. The C++ comparison build uses Linux x86-64, GCC 13 or later and GNU Make:

```sh
make -j3 all
make check
```

The dependency script downloads PEIGEN headers at the commit specified in `third_party/peigen.json` and verifies their hashes.

## Algorithms

| Paper algorithm | Task | Implementation |
|---|---|---|
| 2 | Minimum degree | `Engine::minimum` in `src/optimized/degree.hpp` and `src/minmax/degree.hpp` |
| 3 | Complete degree spectrum | `RankEngine::transform` and `profile` in `src/rank/rank.hpp` |
| 4 | Maximum-degree update | `Engine::apply` with `Metric::Maximum` |
| 5 | Degree-spectrum update | `RankEngine::apply`, `change_coefficients` and `profile` |

All nonzero component vectors are counted, including zero components under the convention degree(0) = 0. Implementation details are in [algorithms.md](docs/algorithms.md).

## Examples

Run [examples/run.sh](examples/run.sh) after building, or follow the individual commands in [examples/README.md](examples/README.md). Examples cover PRESENT, AES and three published S-boxes.

## Experimental data and figures

```sh
make analyze
make figures
```

These commands validate saved results, reconstruct the tables and export the four figures as PDF, SVG and PNG. They do not run timing experiments. Measurement definitions, hardware, workloads and fresh-run commands are in [experiments.md](docs/experiments.md); the constrained search is described in [search.md](docs/search.md).

| Path | Contents |
|---|---|
| `src/` | Algorithms, comparison methods and correctness checks |
| `data/inputs/` | Saved S-box truth tables and transpositions |
| `data/literature/` | Published starting tables, source attribution and resulting tables |
| `data/measurements/minmax/` | 3,520 minimum/maximum timing records |
| `data/measurements/spectrum/` | 2,705 degree-spectrum timings, 75 prefix-cost runs and 24 operation-count records |
| `data/measurements/search/` | 45 search timings, nine warmups, nine diagnostics and nine verification runs |
| `results/` | Reconstructed statistics and verification summaries |
| `figures/` | The four experimental figures |
| `scripts/` | Analysis, plotting, dependency and measurement utilities |
| `examples/` | Runnable usage examples |

The 8-bit comparison of PEIGEN and bitwise recomputation is also available at [results/round2/overlap-8bit.csv](results/round2/overlap-8bit.csv). Download all files through GitHub's **Code → Download ZIP** option.

## License and citation

Project code is GPL-3.0-only; see [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff). PEIGEN retains its own GPL-3.0 license and notices. Please also cite Bao, Guo, Ling and Sasaki, *PEIGEN*, DOI: 10.13154/tosc.v2019.i1.330-394. Published S-box tables retain the attribution and licenses in [data/literature/README.md](data/literature/README.md).
