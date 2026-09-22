#pragma once
#include <cstdint>
#include <ostream>

namespace coam {
// Logical work, not instruction counts. Only the COAM_COUNT build updates these.
#define COAM_COUNTER_FIELDS(X) \
    X(scalar_mobius_xor) X(scalar_component_xor) X(scalar_elimination_xor) \
    X(scalar_coefficient_updates) X(scalar_rollback_xor) X(scalar_degree_checks) \
    X(pivot_checks) X(delta_candidates) X(delta_coefficients) \
    X(word_mobius_xor) X(word_component_xor) X(word_elimination_xor) \
    X(word_coefficient_updates) X(word_rollback_xor) X(word_degree_checks) \
    X(top_guard_checks) X(top_guard_kept) X(top_guard_fallback) X(degree_drops) \
    X(scalar_guard_checks) X(scalar_guard_kept) X(scalar_guard_fallback) X(scalar_degree_drops)
struct Counts {
#define FIELD(name) uint64_t name = 0;
    COAM_COUNTER_FIELDS(FIELD)
#undef FIELD
};
inline Counts counts;
inline void print_counts(std::ostream& out) {
    out << '{';
    bool first = true;
#define FIELD(name) if (!first) out << ','; first = false; out << '"' << #name << "\":" << counts.name;
    COAM_COUNTER_FIELDS(FIELD)
#undef FIELD
    out << '}';
}
}
#ifdef COAM_COUNT
#define COAM_ADD(name, amount) (::coam::counts.name += (amount))
#else
#define COAM_ADD(name, amount) ((void)0)
#endif
