# Experimental protocol

## Outputs and inputs

For an n-by-m mapping, the minimum and maximum degrees are taken over all nonzero component masks. The complete degree histogram is `H[d] = number of nonzero masks whose component has degree d`. Spectrum evaluation returns this complete histogram; the local-search score is its degree sum, `sum(d * H[d])`. All implementations use degree zero for the zero Boolean function. The optional enumeration of masks attaining the minimum is checked for correctness but is outside the scalar minimum timing.

The saved corpus contains ten random permutations at each width 3 through 16, six practical instances, and one 8-bit identity search control. Each LUT begins with `n m`, followed by decimal output values in ascending input order. Swap files contain decimal input-index pairs. The explicit files, rather than a language-dependent random generator, define the workload. Seeds, file digests and source URLs are in [the input manifest](../data/input-manifest.json).

The practical instances are PRESENT (4 bits), AES and CLEFIA S0/S1 (8 bits), MISTY1 S9 (9 bits), and MISTY1 FI with subkey 0 (16 bits). FI is a fixed-key nonlinear component. Initial standard-instance properties are tabulated separately; transposition measurements follow the resulting variants. AES was checked against the fixed PEIGEN LUT; the MISTY1 implementation was checked against both RFC 2994 block test vectors.

## Single-thread suite

Every ordinary timing configuration has five repetitions. Random scaling uses ten independently generated inputs at each width. Repetitions restart from the same saved input in separate processes.

**A: algorithm reference comparison.** Both sides use contiguous `int32_t` ANF coefficients, cached monomial weights, reusable buffers and the same compiler options. `reference` transforms coordinate truth tables and independently constructs components as needed. `reference-proposed` implements the proposed elimination or symmetric-difference update. Minimum timing is repeated static evaluation; maximum/spectrum timing commits every saved transposition. A natural swap profile and a high-update profile are both retained. The latter uses pairs affecting at least half the ANF index set. Calibrated batch sizes were frozen before formal timing and are recorded in `protocol/workloads.json` and the task list. A separate instrumented build replays the same workloads for logical-operation counts; its clocks never enter performance results.

**B: optimized implementation comparison.** Widths 3–8 call the actual PEIGEN numerical kernels through a small adapter. Widths 9–16 use the project's packed traditional implementation (`bitwise`), with 64-bit truth tables, a bitwise Möbius transform and Gray-order component reuse. The proposed implementation uses packed storage, reused buffers, sparse symmetric differences, a highest-degree affected-term check and in-place updates/rollback. Dynamic baselines keep prepared coordinate truth tables and recompute the requested property after each swap. Static timing includes the LUT-to-packed-input conversion at each call. Each update cycle contains 128 fixed transpositions. All cycles are committed in order; the state continues between cycles. The number of cycles depends on width and metric, and is identical across compared methods. For example, the 16-bit B workloads have 1 minimum evaluation, 512 maximum updates or 128 spectrum updates per timed batch.

The packed transform and weight-order degree scan reproduce/adapt Bakoev's [2019 method](https://arxiv.org/abs/1905.08649) and [2020 method](https://arxiv.org/abs/2007.01116). This is a project implementation of those techniques. It uses six within-word transform stages followed by word-wise XOR stages, bitwise weight-lexicographic order for n <= 8 and its CB-WLO variant above 8. Gray-order enumeration reuses one component buffer; the maximum-degree baseline only evaluates coordinates. The PEIGEN adapter uses `get_coordinates_ANF` and `degree_from_ANF` for maximum degree, SIMD `degree` for 3–4-bit minimum/spectrum, and `get_components_ANF` with weight-order `degree_from_ANF` for 5–8-bit minimum/spectrum. Checks compare this adapter against the upstream spectrum interface and exclude the zero mask from the histogram.

**Complete local search.** Fourteen 8-bit starting permutations are used: ten random inputs, AES, CLEFIA S0/S1 and identity. Three separate searches maximize minimum degree, maximum degree or spectrum sum. Candidates are unordered index pairs scanned in lexicographic order. The first strict improvement is accepted and scanning restarts; rejected swaps are undone. Search ends after a complete neighborhood scan with no improvement. Thus a zero-acceptance search still examines 32,640 candidates. The evaluator changes between methods; the start, candidate order, acceptance and stopping rules remain fixed. Saved candidate digests, accepted swaps, final LUTs and complete histograms agree. This is a local-search workload, not an enumeration of every possible S-box.

## Timing and statistics

`context_ns` measures dimension-dependent setup. `init_ns` measures evaluator-state construction, including its first property evaluation. `kernel_ns` measures subsequent static evaluations, committed updates or the complete local search. `total_ns` is exactly their sum. File reading, correctness checks and result writing are outside these intervals. `process_seconds` also includes process overhead, validation and output.

Per-step cost is `kernel_ns / steps`. Batch-total cost includes setup and initialization; it is not a single cold-call latency. In A, LUT conversion and state construction are also recorded separately. Large spectrum initialization can materially reduce total-time speedup.

For single-thread results, take the median of five times for each input, form baseline/proposed ratios within that input, then report the median and interquartile range across ten inputs. Quartiles use linear interpolation over ordered observations. All repetitions and slow configurations are retained. The published scripts recompute all tables and verify saved outputs and work sizes.

## Measurement environment

Measurements used Ubuntu Linux 6.8.0-134-generic, GCC 13.3.0, two AMD EPYC 9754 sockets (256 physical cores / 512 logical CPUs), and approximately 1 TiB of RAM. Build flags were `-std=c++20 -O3 -march=native -fopenmp -Wall -Wextra`. Timed runs bound CPU 0 and NUMA node 0.

The machine was shared. The scalar-reference archive contains 824 group-A commands that overlapped for approximately 44 minutes with a separate workload on the other socket and NUMA node. Optimized group-B, search and counting batches did not overlap that workload. The affected secondary records remain in the archive; the optimized comparison and search results are unaffected.
