# Figure guide and captions

The current set contains **four figures and one numerical table**. Blue identifies the proposed method or its speedup; red identifies the optimized traditional implementation. Bars and speedup panels use zero-based linear axes; execution-time panels use logarithmic axes.

**Timing.** Execution time in Figures 1–3 and Table 1 excludes initialization. Figure 4 reports total execution time, including initialization and the complete search. Speedup is traditional-method time divided by proposed-method time. The archived tables retain both timing measurements.

**Statistics.** Each random input contributes its median over five repetitions. Curves show medians across ten inputs; speedups are computed for each input before taking their median. Bands show the middle 50% across inputs.

Every figure/table has an embedded-font PDF, editable-text SVG and 300 dpi PNG, 7.4 inches wide. Data and export hashes are in [the manifest](../figures/manifest.json).

Rebuild with:

~~~sh
python3 scripts/analyze.py --check
python3 scripts/plot_results.py
~~~

Use `--only 03_degree_spectrum_update` to regenerate one plot, or `--output build/figure-review` to inspect a candidate before replacing published exports.

## 1–3. Optimized single-thread computation

![Minimum degree](../figures/01_minimum_degree.png)

**Figure 1.** Minimum algebraic degree computation on random n-by-n permutations, n = 3–19. Left: execution time per calculation. Right: speedup. Initialization is excluded. PEIGEN is used at 3–8 bits and the bitwise traditional implementation at 9–19 bits; the vertical dotted line marks this change. Ten inputs per size, five repetitions per input; medians and middle-50% bands as defined above.

![Maximum degree update](../figures/02_maximum_degree_update.png)

**Figure 2.** Maximum algebraic degree updates after output transpositions, with the same traditional implementations and statistics as Figure 1. Left: execution time per update. Right: speedup. Initialization is excluded. All transpositions are retained in sequence.

![Algebraic degree spectrum update](../figures/03_degree_spectrum_update.png)

**Figure 3.** Algebraic degree spectrum updates after output transpositions, with the same traditional implementations and statistics as Figure 1. Left: execution time per update. Right: speedup. Initialization is excluded. Each update computes the full algebraic degree spectrum, and all transpositions are retained in sequence.

Source for Figures 1–3: [scaling-results.csv](../results/serial/optimized/scaling-results.csv). Frozen metric-dependent batch sizes are described in [the protocol](protocol.md).

## 4. Published S-box postprocessing at strict NL 104

![Complete search times](../figures/04_search_time.png)

**Figure 4.** Total time to minimum degree 7 for Freyre S8 and Kuznetsov S-box2/S-box3 while strictly preserving NL 104. Each panel has its own time scale. Bars give five-run medians, whiskers give the full range, and dots show all five runs. Total times include initialization, candidate degree evaluation, full-component NL checks and rollback. Both evaluators follow the same trajectory and produce the same final LUT.

Source: [literature.json](../results/literature.json), independently reconstructed from [the 30 raw timing records](../data/measurements/literature.jsonl).

## Table 1. Practical instances

![Practical-instance speedups](../figures/table01_practical_instances.png)

**Table 1.** Speedups on six practical inputs, excluding initialization. Each method contributes its median over five repetitions. Traditional implementations are PEIGEN up to 8 bits and the bitwise method at 9 and 16 bits. FI uses fixed subkey 0. Minimum-degree calculation uses the original input; update measurements follow fixed transpositions starting from it.

Source: [practical-results.csv](../results/serial/optimized/practical-results.csv).
