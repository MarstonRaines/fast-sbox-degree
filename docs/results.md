# Results guide

The archive contains 13,000 serial timing records, 30 literature-postprocessing timings, and 1,460 separate operation-count records. All three degree metrics cover every width from 3 through 19.

## Fixed-workload comparison

| Width | Minimum | Maximum updates | Spectrum updates | Baseline |
|---|---:|---:|---:|---|
| 8 | 10.66× | 12.40× | 10.85× | PEIGEN |
| 16 | 448.29× | 96.07× | 6.79× | Bitwise |
| 17 | 813.65× | 188.65× | 8.83× | Bitwise |
| 18 | 1452.06× | 229.46× | 11.25× | Bitwise |
| 19 | 2728.09× | 319.38× | 15.29× | Bitwise |

These are median paired speedups across ten inputs, using each evaluator’s five-run median and excluding initialization. Each 16–19-bit workload contains one minimum calculation, 512 maximum updates, or 128 full-spectrum updates. The [full table](../results/serial/optimized/scaling-results.csv) also includes initialization and total times. Slower 4-bit cases remain in the archive. The [practical table](../results/serial/optimized/practical-results.csv) covers six named components.

## Published S-boxes at strict NL 104

Three original tables from two papers reach minimum degree 7 and the saturated spectrum of 255 degree-7 components. Every accepted state retains standard vectorial NL 104 and differential uniformity 8.

| Starting table | Degree sum | Candidates | PEIGEN (ms) | Proposed (ms) | Speedup |
|---|---|---:|---:|---:|---:|
| freyre2020-s8 | 1778 → 1785 | 17,555 | 169.234 | 114.472 | 1.478× |
| kuznetsov2023-hill-2 | 1784 → 1785 | 488 | 7.382 | 6.036 | 1.223× |
| kuznetsov2023-hill-3 | 1784 → 1785 | 428 | 5.592 | 4.405 | 1.270× |

Times include initialization, degree checks, shared NL checks, and rollback. Both evaluators use identical candidates, acceptance and stopping rules. Freyre S8 takes three accepted swaps, with degree-7 coefficient ranks 5 → 6 → 7 → 8. A swap changes this block by rank at most one, so this path attains the three-swap lower bound. Strict minimum-only improvement stalls on this starting table; improving the spectrum allows progress before the minimum changes.

The [sources and witnesses](../data/literature/sources-and-witnesses.json) give full original/final LUTs and paths. `python3 scripts/literature.py --check` independently recomputes all component degrees and Walsh transforms, differential uniformity, and ranks.

## Additional archived evidence

Group A and its operation-count tables cover 3–16 bits and support the complexity analysis. The [legacy search table](../results/serial/optimized/search-results.csv) retains 42 unconstrained degree searches from fourteen 8-bit starts using the optimized build. Those runs are separate from the strict-NL experiment shown in Figure 4; they do not establish preservation of NL.
