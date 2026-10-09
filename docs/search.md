# Fixed-NL/DU search

Measurements are in `data/measurements/search/`. The three inputs are published 8-bit S-boxes; each has initial minimum degree 6, nonlinearity 104 and differential uniformity 8.

## Search variants

| Variant | Degree evaluation | Strict improvement objective |
|---|---|---|
| `peigen` | PEIGEN minimum-degree routine | Minimum component degree only |
| `rank-static` | Static rank-based degree spectrum | Lexicographic (minimum degree, spectrum sum) |
| `rank-shared` | Incremental rank-based degree spectrum | The same lexicographic pair |

All use unordered input pairs `(p,q)` in lexicographic order, first strict improvement, restart after an accepted swap, and stop at minimum degree 7 or a full unsuccessful sweep. Every accepted state has NL exactly 104 and differential uniformity exactly 8. The PEIGEN group's score has no spectrum tie-break. These are the search configurations measured here; PEIGEN also provides other evaluation and generation capabilities.

The screening order is objective improvement, NL equality, then DU equality. A candidate already rejected by its score does not need a constraint update; one rejected by NL does not need a DDT update. Each modified cache is restored on rejection. Consequently, the PEIGEN group on Freyre S8 performs no constraint updates: every candidate fails its degree objective, and the initial admissible box remains unchanged.

## Menyachikhin implementation correspondence

Reference: A. V. Menyachikhin, *The change in linear and differential characteristics of substitution after the multiplication by transposition*, Mathematical Aspects of Cryptography 11(2), 111–123 (2020), [DOI 10.4213/mvk325](https://doi.org/10.4213/mvk325).

`src/rank/search_constraints.hpp` implements the paper's Algorithms 1 and 2 (pp. 114 and 117). `menyachikhin_walsh_update` enumerates the input/output masks whose products with the two swap differences are odd. It changes each affected unnormalised Walsh entry by signed 4 and updates its absolute-value frequency. `menyachikhin_ddt_update` skips the input difference `p xor q`; equal old output differences use the two-entry ±4 branch, otherwise the four-entry ±2 branch. Both maintain exact maxima through frequencies.

The paper's LAT and DDT use normalised values. Our integer tables are multiplied by `2^n`; NL is `2^(n-1) - max(abs(W))/2`, and DU is the maximum DDT count. Invariant zero-mask entries are stored too; they do not affect either maximum for a bijection. No full Walsh transform or DDT recomputation runs inside candidate checks. Initial tables and independent validation use full calculations.

The published Examples 1 and 2 are explicit tests in `src/rank/check_constraints.cpp`. The same tests cover 3,154 swaps, continuous commits, rejected-candidate rollback, same-index swaps and changes in both directions of each maximum. ASan/UBSan and checked container bounds pass. Nine verification searches pair every candidate degree with PEIGEN's complete spectrum, recompute every changed Walsh/DDT table and verify rollback. Independent ANF, direct Walsh sums and DDT checks verify every accepted state outside timing.

## Measurements and results

Measurements use a shared dual AMD EPYC 9754 machine with GCC 13.3.0. The hardware and compiler settings are recorded in `data/environment.json`; the raw records retain measured binary and input hashes. Timing uses one thread, CPU 0, NUMA node 0, a 50 ms arithmetic warmup and one whole-search warmup per input/group. Five fresh processes per group rotate execution order. There are 45 formal records, nine verification searches, nine warmups and nine diagnostic searches. No claim of whole-machine exclusivity is made.

`total_ns = context_ns + init_ns + kernel_ns`, with `init_ns = degree_init_ns + constraints_init_ns`. Total time includes both constraint initialisation and every actual constraint update and rollback. Input reading, CPU warmup, result serialization and independent correctness checks are outside timing. Separate diagnostic clocks do not contaminate the formal records.

| Input | PEIGEN | Static spectrum | Spectrum update |
|---|---|---|---|
| Freyre S8 | 104.061 ms; stops at 6; 32,640 candidates | 296.029 ms; reaches 7; 17,555 candidates | 287.126 ms; reaches 7; 17,555 candidates |
| Kuznetsov S-box2 | 9.949 ms; reaches 7; 488 candidates | 8.495 ms; reaches 7; 488 candidates | 8.316 ms; reaches 7; 488 candidates |
| Kuznetsov S-box3 | 9.055 ms; reaches 7; 428 candidates | 7.808 ms; reaches 7; 428 candidates | 7.694 ms; reaches 7; 428 candidates |

Freyre's stopping time is not a time to degree 7, and no PEIGEN/update speedup is reported for that input. Static/update paths agree exactly. On the two Kuznetsov inputs the PEIGEN group has the same accepted path as the spectrum groups. The eight distinct initial/accepted states independently recompute to NL 104 and DU 8. DU rejects no additional candidate in these three examples after the objective and NL filters; its equality condition is nevertheless explicitly evaluated before every accepted swap.

## Reproduction

Regenerate statistics and Figure 4 without rerunning experiments:

```sh
python3 scripts/analyze_search.py
python3 scripts/plot_search.py
```

Statistics and independently recomputed states are under `results/search/`; Figure 4 is under `figures/`. Use `python3 scripts/run_measurements.py search` for the saved search protocol. Ordinary correctness checks use `make rank-check`.
