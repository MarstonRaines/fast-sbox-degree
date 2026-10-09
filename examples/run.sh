#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
build/minmax-optimized proposed min data/inputs/PRESENT.lut data/inputs/PRESENT.swaps 16
build/minmax-optimized proposed max data/inputs/AES.lut data/inputs/AES.swaps 16
build/rank-optimized bench rank-static data/inputs/AES.lut data/inputs/AES.swaps 16
build/rank-optimized bench rank-shared data/inputs/AES.lut data/inputs/AES.swaps 16
build/rank-search verify rank-shared data/literature/inputs/kuznetsov2023-hill-2.lut
