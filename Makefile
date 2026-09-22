CXX = g++
PYTHON ?= python3
CXXFLAGS ?= -std=c++20 -O3 -march=native -fopenmp -Wall -Wextra
CPPFLAGS += -Ithird_party
OPTFLAGS = -DCOAM_OUTPUT_PACKED_MINIMUM -DCOAM_TILED_SPECTRUM -DCOAM_SPARSE_DELTA

.PHONY: all dependencies check check-full analyze figures clean
all: build/serial build/counts build/verify-full build/optimized build/extended build/minmax build/postprocess build/warm_cpu build/verify-minmax

dependencies: third_party/PEIGEN/.verified
third_party/PEIGEN/.verified: third_party/peigen.json scripts/fetch_peigen.py
	$(PYTHON) scripts/fetch_peigen.py

build:
	mkdir -p build

build/serial: $(wildcard src/serial/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/serial/bench.cpp src/serial/peigen_adapter.cpp -o $@

build/counts: $(wildcard src/serial/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOAM_COUNT src/serial/bench.cpp src/serial/peigen_adapter.cpp -o $@

build/optimized build/extended build/minmax: build/%: src/%/bench.cpp src/%/degree.hpp src/%/reference.hpp src/%/counts.hpp src/%/peigen_adapter.cpp third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) src/$*/bench.cpp src/$*/peigen_adapter.cpp -o $@

build/postprocess: src/postprocess.cpp $(wildcard src/optimized/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/optimized src/postprocess.cpp src/optimized/peigen_adapter.cpp -o $@

build/warm_cpu: src/warm_cpu.cpp | build
	$(CXX) -std=c++20 -O3 $< -o $@

build/verify-minmax: src/verify_minmax.cpp $(wildcard src/minmax/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -Isrc/minmax src/verify_minmax.cpp src/minmax/peigen_adapter.cpp -o $@

build/verify-full: src/verify_full_state.cpp $(wildcard src/serial/*.hpp) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc/serial src/verify_full_state.cpp src/serial/peigen_adapter.cpp -o $@

check: all
	OMP_NUM_THREADS=1 build/serial check
	OMP_NUM_THREADS=1 build/counts check
	OMP_NUM_THREADS=1 build/optimized check
	OMP_NUM_THREADS=1 build/verify-minmax --edges
	$(PYTHON) scripts/analyze.py --check
	$(PYTHON) scripts/literature.py --check

check-full: build/verify-full
	build/verify-full data/inputs/random_n16_00.lut data/inputs/MISTY1_FI_key0000.lut

analyze:
	$(PYTHON) scripts/analyze.py

figures:
	$(PYTHON) scripts/plot_results.py

clean:
	rm -rf build
