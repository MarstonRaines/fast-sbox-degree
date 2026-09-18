# Results guide

All speedups are baseline time divided by proposed-method time, unless explicitly labeled hardware scaling. Values below summarize the corresponding complete CSV tables; they do not replace the per-input records.

## Optimized single-thread comparison

| Width / baseline | Minimum evaluation | Maximum update | Spectrum update |
|---|---:|---:|---:|
| 8 / PEIGEN | 1.21x | 8.66x | 3.04x |
| 9 / packed traditional | 1.35x | 5.37x | 2.48x |
| 16 / packed traditional | 25.91x | 37.14x | 2.52x |

These are median per-input core-time ratios across ten random permutations, with five repeats per input. At 16 bits, the corresponding batch-total ratios are 18.64x, 15.50x and 2.05x. Work sizes are 1 static evaluation, 512 maximum updates and 128 spectrum updates. Some small-input cases favor PEIGEN; they remain in [the full scaling table](../results/serial/optimized/scaling-results.csv). [Practical-instance results](../results/serial/optimized/practical-results.csv) keep the six named cases separate from random-input summaries.

## Reference comparison and operation counts

At 16 bits, minimum evaluation is 5345.13x faster than the specified scalar reference. High-update maximum and spectrum workloads give 46.51x and 48.88x core-time ratios. Natural swaps are usually sparser and give larger ratios. For the high-update spectrum workload, the logical ANF XOR count ratio is 32.0078125, close to the reference growth factor 2n = 32. Logical operation counts and asymptotic factors are not exact runtime predictions: scans, branches, memory access and initialization also contribute.

Use [reference-scaling.csv](../results/serial/reference-scaling.csv) together with [mechanism-per-input.csv](../results/serial/mechanism-per-input.csv). These use the scalar representation described in the protocol. Their ratios are separate from the PEIGEN/packed comparison.

## Complete 8-bit search

For ten random starting permutations, median per-input total-time ratios against PEIGEN are **7.03x / 6.26x / 2.06x** for minimum degree, maximum degree and spectrum sum. Both implementations follow the same search trajectory and obtain the same final result. Five minimum-degree searches improve their objective; seven spectrum searches improve theirs. Random maximum-degree starts already have degree 7 and still complete the full neighborhood scan.

For CLEFIA S0, spectrum sum increases from 1530 to 1785, with 8 accepted swaps and 32,899 examined candidates. Median total times are 112.86 ms for PEIGEN and 56.25 ms for the proposed evaluator. All fourteen starts, their initial/final scores, candidate counts and times are in [search-results.csv](../results/serial/optimized/search-results.csv); complete accepted swaps and final permutations are linked from the raw records.

## Parallel extension

For the two fixed random 16-bit inputs, the proposed minimum evaluator achieves **4.13x / 4.18x core hardware speedup at 16 threads**, with **2.13x / 2.15x batch-total speedup**. Spectrum updates achieve **17.32x / 16.69x core speedup at 64 threads**, with **4.59x / 4.71x batch-total speedup**. MISTY1 FI spectrum updates give 18.61x core and 9.22x total speedup at 64 threads.

At the same 64-thread cap, spectrum core-time ratios of the packed traditional baseline to the proposed method are 8.66x / 7.82x for the two random inputs; total-time ratios are 1.86x / 1.81x. Small 8-bit search workloads obtain approximately 1x hardware speedup because their kernels mostly stay below the parallel threshold. Increasing the cap does not monotonically improve every task.

Use [timings-and-scaling.csv](../results/parallel/timings-and-scaling.csv) for all absolute times, [hardware-scaling.csv](../results/parallel/hardware-scaling.csv) for paired core/total ratios and quartiles, [algorithm-comparison.csv](../results/parallel/algorithm-comparison.csv) for equal-cap comparisons, and [serial-anchors.csv](../results/parallel/serial-anchors.csv) for measured serial-build differences. Parallel A uses one random input per width and B uses two; this is a scaling subset, distinct from the ordinary ten-input study.
