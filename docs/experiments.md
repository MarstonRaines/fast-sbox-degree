# Experimental measurements

The measurements use Ubuntu Linux, GCC 13.3.0, C++20, `-O3 -march=native -fopenmp`, one thread, CPU 0 and NUMA node 0 on a shared dual AMD EPYC 9754 system. Each invocation has an untimed 50 ms arithmetic warmup. Memory is Linux VmHWM after exec; MiB denotes 2^20 bytes. `data/environment.json` records the settings, and each measurement retains input and binary hashes.

## Degree calculations

The saved suite contains 170 random bijections (ten per size from 3 through 19 bits) and six practical inputs, including PRESENT, AES, CLEFIA S0/S1 and MISTY1 components. Each input has five repetitions per method. PEIGEN is used at 3–8 bits, and bitwise recomputation at 9–19 bits; both recomputation baselines are also evaluated on the 13 eight-bit inputs.

There are 3,520 minimum/maximum measurements and 2,705 complete-spectrum measurements. `protocol.json` in each data directory gives every input, method, workload and repetition. Every static evaluation includes the truth-table-to-ANF transform. Update kernels exclude one-time context construction and evaluator initialization; those costs are recorded separately. All methods use the same saved transposition stream: a prefix for fewer than 128 swaps, cyclic continuation otherwise. The spectrum implementations compute the same complete multiplicities for every compared state.

Statistics first take the median of five repetitions within each input/method. Speedups compare methods on the same input; the summary then takes the median and interquartile range over the ten random inputs at each size. Practical cases are listed individually. No individual timing is removed from these datasets.

## Memory, cumulative cost and operation counts

The spectrum records contain context/state allocation and whole-process peak resident memory. The 75 prefix runs measure cumulative costs, including preprocessing and initialization, at K = 1, 16, 128 and 512; 19-bit runs end at 128. Each method uses the same input and checkpoints. The 24 untimed diagnostic runs count affected coefficients, scanned rows and basis XORs at K = 0, 1, 16 and 128. Checkpoint instrumentation is confined to the prefix runs.

## Search

The three starting tables and attribution are in `data/literature/`. All accepted states have nonlinearity 104 and differential uniformity 8. PEIGEN scores minimum degree only; the two spectrum methods score the lexicographic pair (minimum degree, spectrum sum). Total time includes preprocessing, degree and constraint initialization, candidate evaluation, constraint updates and rollback. Input reading and independent verification are outside timing. See [search.md](search.md) for the complete protocol and results, including Freyre S8 where the PEIGEN group stops at degree 6.

## Reconstruct the saved results

```sh
make analyze
make figures
```

The analyzers check every recorded repetition, input hash and paired final output. Search analysis independently recomputes the ANF, Walsh transform and difference distribution table of each distinct initial/accepted state. `data/manifest.json` provides checksums for the data files. Result CSVs use `results/minmax/`, `results/spectrum/` and `results/search/`.

## Run measurements

Build with `make -j3 all diagnostics`. Inspect the command queue with:

```sh
python3 scripts/run_measurements.py minmax --dry-run
python3 scripts/run_measurements.py spectrum --dry-run
python3 scripts/run_measurements.py supplemental --dry-run
python3 scripts/run_measurements.py search --dry-run
```

Remove `--dry-run` to execute a batch. The runner uses the saved protocol jobs and writes into `runs/<batch>/`; it refuses to overwrite that directory. CPU/NUMA binding defaults to 0 and can be selected with `--cpu` and `--numa`. Run only one timing batch at a time. A different machine or build produces its own measurements and binary hashes.
