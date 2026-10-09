# Usage examples

Build from the repository root with `make -j3 all`. Then run `sh examples/run.sh`, or use the commands below individually.

## Minimum and maximum degree

```sh
build/minmax-optimized proposed min data/inputs/PRESENT.lut data/inputs/PRESENT.swaps 16
build/minmax-optimized proposed max data/inputs/AES.lut data/inputs/AES.swaps 16
```

The minimum-degree command evaluates the saved truth table. The maximum-degree command commits the specified number of saved transpositions. Output is JSON, including the degree, elapsed times and output digest.

## Degree spectrum

```sh
build/rank-optimized bench rank-static data/inputs/AES.lut data/inputs/AES.swaps 16
build/rank-optimized bench rank-shared data/inputs/AES.lut data/inputs/AES.swaps 16
```

Both commands return the same final histogram and truth-table digest. `rank-static` recomputes the spectrum for each candidate; `rank-shared` updates it after each transposition.

## Constrained S-box search

```sh
build/rank-search verify rank-shared data/literature/inputs/kuznetsov2023-hill-2.lut
```

This checks all candidate spectra and the NL/DU constraint updates. The output includes the accepted transpositions and final truth table. Use `peigen` or `rank-static` for the other search configurations described in [search.md](../docs/search.md).
