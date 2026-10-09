CXX = g++
PYTHON ?= python3
CXXFLAGS ?= -std=c++20 -O3 -march=native -fopenmp -Wall -Wextra
CPPFLAGS += -Ithird_party
OPTFLAGS = -DCOAM_OUTPUT_PACKED_MINIMUM -DCOAM_SPARSE_DELTA
.PHONY: all dependencies check analyze figures rank-all rank-check
all: rank-all build/minmax-optimized build/minmax-minmax build/verify-minmax build/verify-minmax-optimized build/verify-small build/verify-small5
dependencies: third_party/PEIGEN/.verified
third_party/PEIGEN/.verified: third_party/peigen.json scripts/fetch_peigen.py
	$(PYTHON) scripts/fetch_peigen.py
build:
	mkdir -p build
rank-all: build/rank-optimized build/rank-extended build/rank-search
build/rank-optimized build/rank-extended: build/rank-%: $(wildcard src/rank/*.hpp) src/rank/bench.cpp src/%/degree.hpp src/%/counts.hpp src/%/peigen_adapter.cpp third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/rank -Isrc/$* src/rank/bench.cpp src/$*/peigen_adapter.cpp -o $@
build/rank-search: $(wildcard src/rank/*) $(wildcard src/optimized/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/rank -Isrc/optimized src/rank/search.cpp src/optimized/peigen_adapter.cpp -o $@
build/minmax-optimized build/minmax-minmax: build/minmax-%: src/minmax_bench.cpp src/%/degree.hpp $(wildcard src/rank/*.hpp) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -DCOAM_SMALL_MINIMUM_5 -Isrc/rank -Isrc/$* src/minmax_bench.cpp src/$*/peigen_adapter.cpp -o $@
build/verify-minmax: src/verify_minmax.cpp src/minmax/degree.hpp $(wildcard src/rank/*.hpp) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -DCOAM_SMALL_MINIMUM_5 -Isrc/rank -Isrc/minmax src/verify_minmax.cpp src/minmax/peigen_adapter.cpp -o $@
build/verify-minmax-optimized: src/verify_minmax.cpp src/optimized/degree.hpp src/optimized/peigen_adapter.cpp $(wildcard src/rank/*.hpp) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -DCOAM_SMALL_MINIMUM_5 -Isrc/rank -Isrc/optimized src/verify_minmax.cpp src/optimized/peigen_adapter.cpp -o $@
build/adapter-optimized.o build/adapter-extended.o: build/adapter-%.o: src/%/peigen_adapter.cpp src/%/degree.hpp third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/$* -c $< -o $@
build/verify-small5: src/verify_small5.cpp src/rank/small_min5.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
build/verify-small: src/verify_small.cpp src/rank/small_rank.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@
build/check-constraints: src/rank/check_constraints.cpp src/rank/search_constraints.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc/rank -Isrc/optimized $< -o $@
rank-check: rank-all build/check-constraints
	build/check-constraints data/literature/inputs/freyre2020-s8.lut data/literature/inputs/kuznetsov2023-hill-2.lut data/literature/inputs/kuznetsov2023-hill-3.lut
	build/rank-search verify peigen data/literature/inputs/freyre2020-s8.lut
	build/rank-search verify peigen data/literature/inputs/kuznetsov2023-hill-2.lut
	build/rank-search verify peigen data/literature/inputs/kuznetsov2023-hill-3.lut
	build/rank-optimized check
	build/rank-search verify rank-shared data/literature/inputs/freyre2020-s8.lut
	build/rank-search verify rank-shared data/literature/inputs/kuznetsov2023-hill-2.lut
	build/rank-search verify rank-shared data/literature/inputs/kuznetsov2023-hill-3.lut
	build/rank-search verify rank-static data/literature/inputs/freyre2020-s8.lut
	build/rank-search verify rank-static data/literature/inputs/kuznetsov2023-hill-2.lut
	build/rank-search verify rank-static data/literature/inputs/kuznetsov2023-hill-3.lut
check: all
	build/verify-minmax
	build/verify-minmax-optimized
	build/verify-small
	build/verify-small5
	$(MAKE) rank-check
analyze:
	$(PYTHON) scripts/analyze_minmax.py
	$(PYTHON) scripts/analyze_spectrum.py
	$(PYTHON) scripts/analyze_search.py
figures:
	$(PYTHON) scripts/plot_minmax.py
	$(PYTHON) scripts/plot_spectrum.py
	$(PYTHON) scripts/plot_search.py

.PHONY: diagnostics
diagnostics: build/prefix-optimized build/prefix-extended build/rank-counts
build/prefix-optimized build/prefix-extended: build/prefix-%: src/diagnostics/prefix.cpp build/adapter-%.o $(wildcard src/rank/*.hpp) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/rank -Isrc/$* $< build/adapter-$*.o -o $@
build/rank-counts: src/rank/bench.cpp build/adapter-extended.o $(wildcard src/rank/*.hpp) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -DCOAM_COUNT -Isrc/rank -Isrc/extended $< build/adapter-extended.o -o $@
