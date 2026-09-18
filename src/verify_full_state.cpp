// Additional untimed check: full packed coefficient state, not just degrees.
#include "degree.hpp"
#include <fstream>
#include <iostream>
#include <memory>

using namespace coam;

Box read_box(const char* path, unsigned& n, unsigned& m) {
    std::ifstream input(path);
    if (!(input >> n >> m)) throw std::runtime_error("input header");
    if (n!=16 || m!=16) throw std::runtime_error("expected 16x16 inputs");
    Box box(1U << n);
    for (auto& value:box) if (!(input>>value)) throw std::runtime_error("input entry");
    return box;
}

void check_full(const Engine& e) {
    // Linearity lets us compare every stored ANF word against an independently built
    // XOR of coordinate ANFs, without reusing the parent recurrence.
    std::vector<Word> coordinates(e.c.m*e.c.words);
    e.c.coordinates(e.box,coordinates);
    std::vector<Word> combination(e.c.words);
    Hist histogram{};
    for (unsigned lambda=1; lambda<e.c.components; ++lambda) {
        std::fill(combination.begin(),combination.end(),0);
        for (unsigned bits=lambda; bits; bits&=bits-1) {
            const Word* row=coordinates.data()+std::countr_zero(bits)*e.c.words;
            for (size_t w=0; w<e.c.words; ++w) combination[w]^=row[w];
        }
        const Word* stored=e.component_anfs.data()+size_t(lambda)*e.c.words;
        if (!std::equal(combination.begin(),combination.end(),stored))
            throw std::runtime_error("full component ANF mismatch at lambda="+std::to_string(lambda));
        unsigned degree=e.c.degree(combination.data());
        ++histogram[degree];
        if (e.degrees[lambda]!=degree)
            throw std::runtime_error("component degree mismatch");
    }
    if (histogram!=e.hist || spectrum_sum(histogram)!=e.value)
        throw std::runtime_error("full state histogram/value mismatch");
}

int main(int argc,char** argv) {
    try {
        if (argc!=3) throw std::runtime_error("requires two saved 16-bit LUTs");
        omp_set_dynamic(0);
        unsigned states=0;
        for (int input=1; input<argc; ++input) {
            unsigned n,m;
            Box initial=read_box(argv[input],n,m);
            Context context(n,m);
            for (int threads:{1,4,64}) {
                omp_set_num_threads(threads);
                Engine engine(context,Method::Proposed,Metric::Spectrum,initial);
                omp_set_num_threads(1);
                check_full(engine); ++states;
                for (auto [p,q]:{std::pair{0U,1U},std::pair{37U,179U}}) {
                    Box before=engine.box;
                    Hist histogram=engine.hist;
                    omp_set_num_threads(threads);
                    engine.apply(p,q,true);
                    omp_set_num_threads(1);
                    check_full(engine); ++states;
                    omp_set_num_threads(threads);
                    engine.undo(p,q);
                    if (engine.box!=before || engine.hist!=histogram)
                        throw std::runtime_error("rollback mismatch");
                    omp_set_num_threads(1);
                    check_full(engine); ++states;
                    omp_set_num_threads(threads);
                    engine.apply(p,q,false);
                }
                omp_set_num_threads(1);
                check_full(engine); ++states;
            }
        }
        std::cout << "{\"status\":\"pass\",\"threads\":[1,4,64],\"inputs\":2,\"full_states\":"
                  << states << ",\"coefficients_per_state\":4294901760,\"includes_update_and_rollback\":true}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
