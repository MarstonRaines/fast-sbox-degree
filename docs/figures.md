# Figure guide and captions

The current set contains **four figures and one numerical table**. Blue identifies the proposed method or its core speedup; red identifies the optimized baseline; gray identifies the speedup including initialization. Direct annotations give measured ratios. Ratios below 1 are retained. Bars and speedup panels use zero-based linear axes; computation-time panels use logarithmic axes. No axis breaks or fitted curves are used.

**Timing.** Core time excludes context and state initialization. "Including initialization" is the ratio of complete fixed-batch times: context + state initialization + timed computation. It is not a single cold-call latency. Speedup is baseline time divided by proposed-method time.

**Statistics.** Each random input first contributes its median over five repetitions. Curves show medians across ten inputs; speedups are computed within each input before summarizing. Bands show the middle 50% (IQR) across those inputs. An annotation is a rounded value at a selected point, not an additional measurement.

Every figure/table has an embedded-font PDF, editable-text SVG and 300 dpi PNG, 7.4 inches wide. Data and export hashes are in [the manifest](../figures/manifest.json).

Rebuild with:

~~~sh
python3 scripts/analyze.py --check
python3 scripts/plot_results.py
~~~

Use `--only 03_degree_spectrum_update` to regenerate one plot, or `--output build/figure-review` to inspect a candidate before replacing published exports.

## 1–3. Optimized single-thread computation

![Minimum degree](../figures/01_minimum_degree.png)

**Figure 1.** Minimum algebraic degree computation on random n-by-n permutations, n = 3–16. Left: core time per evaluation, lower is better. Right: speedup over the optimized baseline, higher is better; gray includes initialization of the full timed batch. PEIGEN is used at 3–8 bits and the packed bitwise traditional implementation at 9–16 bits; the vertical dotted line marks this baseline switch. Ten inputs per size, five repeats per input; median and IQR as defined above. At 16 bits, median core and batch-total speedups are 25.91× and 18.64×.

![Maximum degree update](../figures/02_maximum_degree_update.png)

**Figure 2.** Maximum algebraic degree updates after committed output transpositions, with the same baseline selection and statistics as Figure 1. Left: core time per update. Right: core and initialization-inclusive batch speedups. The 16-bit batch contains 512 updates; its core and batch-total ratios are 37.14× and 15.50×.

![Complete degree-spectrum update](../figures/03_degree_spectrum_update.png)

**Figure 3.** Complete algebraic degree-spectrum updates after committed output transpositions, with the same baseline selection and statistics as Figure 1. Each update maintains the complete component-degree histogram. The 16-bit batch contains 128 updates; its core and batch-total ratios are 2.52× and 2.05×. The full scale-dependent variation, including the 4-bit result below 1 and the 14-bit decrease, is retained.

Source for Figures 1–3: [scaling-results.csv](../results/serial/optimized/scaling-results.csv). Frozen metric-dependent batch sizes are described in [the protocol](protocol.md).

## 4. Complete 8-bit search

![Complete search times](../figures/04_search_time.png)

**Figure 4.** Total time, including initialization, to complete deterministic first-improvement local searches from ten random 8-bit starting permutations. Minimum degree, maximum degree and spectrum sum are optimized in separate searches. Bar heights are medians across per-input five-repeat medians; whiskers show their IQR and dots show all ten per-input medians. Direct time labels give the bar heights. The annotations 7.03×, 6.26× and 2.06× are medians of paired per-input speedups, not ratios recomputed from the displayed bar heights. Both implementations examine the same candidates, accept the same transpositions and return the same final LUT and histogram. Zero-acceptance runs still complete the entire neighborhood scan. Practical and identity-control starts remain in the full table.

Source: [search-results.csv](../results/serial/optimized/search-results.csv).

## Table 1. Practical instances

![Practical-instance speedups](../figures/table01_practical_instances.png)

**Table 1.** Optimized-baseline/proposed speedups on six practical inputs. Each cell gives core speedup followed by initialization-inclusive batch speedup in parentheses. Each method contributes its median over five repeats. Baselines are PEIGEN up to 8 bits and the packed bitwise traditional implementation at 9 and 16 bits. FI is evaluated at fixed subkey 0. Minimum-degree evaluation uses the original input; update measurements follow fixed transpositions starting from it. Ratios below 1 favor the baseline. The approximately 992× minimum-degree ratio on the structured FI instance is specific to that instance; random 16-bit inputs have the separate 25.91× median in Figure 1.

Source: [practical-results.csv](../results/serial/optimized/practical-results.csv).
