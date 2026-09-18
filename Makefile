CXX = g++
PYTHON ?= python3
CXXFLAGS ?= -std=c++20 -O3 -march=native -fopenmp -Wall -Wextra
CPPFLAGS += -Ithird_party

.PHONY: all dependencies check check-full analyze figures clean
all: build/serial build/parallel build/counts build/verify-full

dependencies: third_party/PEIGEN/.verified
third_party/PEIGEN/.verified: third_party/peigen.json scripts/fetch_peigen.py
	$(PYTHON) scripts/fetch_peigen.py

build:
	mkdir -p build

build/serial: $(wildcard src/serial/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/serial/bench.cpp src/serial/peigen_adapter.cpp -o $@

build/parallel: $(wildcard src/parallel/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/parallel/bench.cpp src/parallel/peigen_adapter.cpp -o $@

build/counts: $(wildcard src/serial/*) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOAM_COUNT src/serial/bench.cpp src/serial/peigen_adapter.cpp -o $@

build/verify-full: src/verify_full_state.cpp $(wildcard src/parallel/*.hpp) third_party/PEIGEN/.verified | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Isrc/parallel src/verify_full_state.cpp src/parallel/peigen_adapter.cpp -o $@

check: all
	OMP_NUM_THREADS=1 build/serial check
	OMP_NUM_THREADS=1 build/counts check
	OMP_NUM_THREADS=4 COAM_OMP_MIN_WORK=0 build/parallel check
	$(PYTHON) scripts/analyze.py --check

check-full: build/verify-full
	OMP_WAIT_POLICY=PASSIVE build/verify-full data/inputs/random_n16_00.lut data/inputs/MISTY1_FI_key0000.lut

analyze:
	$(PYTHON) scripts/analyze.py

figures:
	$(PYTHON) scripts/plot_results.py

clean:
	rm -rf build
