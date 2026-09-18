#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <omp.h>

#ifdef COAM_COUNT
#error "experiment3 does not use the serial, shared operation counters"
#endif

namespace coam {
// Small jobs stay serial; count is a ceiling, never an SMT CPU count.
inline uint64_t parallel_threshold() {
    static const uint64_t threshold = [] {
        const char* value = std::getenv("COAM_OMP_MIN_WORK");
        return value ? std::strtoull(value, nullptr, 10) : 32768ULL;
    }();
    return threshold;
}
inline unsigned team_size(unsigned count, uint64_t work) {
    if (work < parallel_threshold()) return 1;
    return std::max(1U, std::min(count, unsigned(omp_get_max_threads())));
}
template<class F> void parallel_for(unsigned first, unsigned last, uint64_t work, F body) {
    unsigned team = team_size(last-first, work);
    if (team == 1) {
        for (unsigned i = first; i < last; ++i) body(i);
    } else {
        #pragma omp parallel for num_threads(team) schedule(static)
        for (unsigned i = first; i < last; ++i) body(i);
    }
}
// Each numbered block owns its scratch and histogram, independently of which
// OpenMP worker executes it. A smaller runtime team cannot omit any block.
template<class F> unsigned parallel_blocks(unsigned first, unsigned last, uint64_t work, F body) {
    unsigned blocks = team_size(last-first, work);
    parallel_for(0, blocks, work, [&](unsigned block) {
        unsigned begin = first + uint64_t(last-first)*block/blocks;
        unsigned end = first + uint64_t(last-first)*(block+1)/blocks;
        body(begin, end, block);
    });
    return blocks;
}
struct alignas(64) Partial {
    std::array<int64_t,17> histogram{};
    unsigned minimum = 0;
};
} // namespace coam
