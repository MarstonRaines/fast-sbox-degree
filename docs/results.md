# Results guide

Speedup is traditional-method execution time divided by proposed-method execution time. The tables below summarize the complete CSV records.

## Optimized single-thread comparison

| Width / baseline | Minimum evaluation | Maximum update | Spectrum update |
|---|---:|---:|---:|
| 8 / PEIGEN | 1.21x | 8.66x | 3.04x |
| 9 / packed traditional | 1.35x | 5.37x | 2.48x |
| 16 / packed traditional | 25.91x | 37.14x | 2.52x |

These speedups exclude initialization. For each input, the execution time is the median of five repetitions; the displayed value is the median speedup across ten random inputs. The [full scaling table](../results/serial/optimized/scaling-results.csv) also reports total execution time, including initialization. The 16-bit workloads contain 1 minimum-degree calculation, 512 maximum-degree updates or 128 spectrum updates. [Practical-instance results](../results/serial/optimized/practical-results.csv) report the six named inputs separately.

## Reference audit tables

The scalar-reference measurements and logical operation counts are retained as secondary audit material for Propositions 2--4. They use the scalar representation described in the protocol and are separate from the practical PEIGEN/packed comparison above. See [reference-scaling.csv](../results/serial/reference-scaling.csv) and [mechanism-per-input.csv](../results/serial/mechanism-per-input.csv).

## Complete 8-bit search

For ten random starting permutations, median per-input total-time ratios against PEIGEN are **7.03x / 6.26x / 2.06x** for minimum degree, maximum degree and spectrum sum. Both implementations follow the same search trajectory and obtain the same final result. Five minimum-degree searches improve their objective; seven spectrum searches improve theirs. Random maximum-degree starts already have degree 7 and still complete the full neighborhood scan.

For CLEFIA S0, spectrum sum increases from 1530 to 1785, with 8 accepted swaps and 32,899 examined candidates. Median total times are 112.86 ms for PEIGEN and 56.25 ms for the proposed evaluator. All fourteen starts, their initial/final scores, candidate counts and times are in [search-results.csv](../results/serial/optimized/search-results.csv); complete accepted swaps and final permutations are linked from the raw records.
