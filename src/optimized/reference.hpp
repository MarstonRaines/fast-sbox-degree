#pragma once
#include "degree.hpp"
#include <chrono>

namespace coam {
// Group A: Algorithm 1's independent combinations versus Algorithms 2/3/4.
// Both sides use int32 coefficients, precomputed weights and reusable buffers.
struct ReferenceEngine {
    const Context& c;
    Method method;
    Metric metric;
    Box box;
    std::vector<int32_t> truth, coord, work, combination, component_anfs;
    std::vector<uint8_t> degrees, old_degrees;
    std::vector<unsigned> delta, top_delta;
    Hist hist{}, old_hist{};
    uint64_t value = 0, old_value = 0, conversion_ns = 0;
    unsigned delta_output = 0;

    ReferenceEngine(const Context& cx, Method me, Metric mt, const Box& s, bool = false)
        : c(cx), method(me), metric(mt), box(s), truth(c.m*c.length),
          coord(c.m*c.length), work(c.m*c.length), combination(c.length) {
        if (method != Method::Reference && method != Method::ReferenceProposed)
            throw std::runtime_error("reference engine requires a group A method");
        if (box.size() != c.length || std::any_of(box.begin(),box.end(),
            [&](uint32_t x) { return x >= c.components; }))
            throw std::runtime_error("invalid truth table");
        auto begin = std::chrono::steady_clock::now();
        for (unsigned b = 0; b < c.m; ++b)
            for (unsigned x = 0; x < c.length; ++x)
                truth[b*c.length+x] = (box[x] >> b) & 1U;
        conversion_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now()-begin).count();
        delta.reserve(c.length);
        top_delta.reserve(c.n);
        if (method == Method::ReferenceProposed && metric != Metric::Minimum) {
            make_coordinates();
            if (metric == Metric::Maximum) {
                degrees.resize(c.m);
                for (unsigned b = 0; b < c.m; ++b) {
                    degrees[b] = ordered_degree(coord.data()+b*c.length);
                    value = std::max(value,uint64_t(degrees[b]));
                }
            } else {
                component_anfs.resize(size_t(c.components)*c.length);
                degrees.resize(c.components);
                for (unsigned l = 1; l < c.components; ++l) {
                    unsigned bit = std::countr_zero(l), parent = l & (l-1);
                    auto* dst = component_anfs.data()+size_t(l)*c.length;
                    const auto* src = component_anfs.data()+size_t(parent)*c.length;
                    const auto* row = coord.data()+bit*c.length;
                    for (unsigned u = 0; u < c.length; ++u) dst[u] = src[u] ^ row[u];
                    COAM_ADD(scalar_component_xor,c.length);
                    degrees[l] = ordered_degree(dst);
                    ++hist[degrees[l]];
                }
                value = spectrum_sum(hist);
            }
            old_degrees.resize(degrees.size());
        } else value = evaluate();
    }

    void make_coordinates() {
        std::copy(truth.begin(),truth.end(),coord.begin());
        for (unsigned b = 0; b < c.m; ++b) {
            auto* row = coord.data()+b*c.length;
            for (unsigned step = 1; step < c.length; step *= 2)
                for (unsigned base = 0; base < c.length; base += 2*step)
                    for (unsigned j = 0; j < step; ++j) row[base+step+j] ^= row[base+j];
        }
        COAM_ADD(scalar_mobius_xor,uint64_t(c.m)*c.n*(c.length/2));
    }

    unsigned full_degree(const int32_t* row) const {
        unsigned d = 0;
        for (unsigned u = 0; u < c.length; ++u)
            if (row[u]) d = std::max(d,unsigned(c.weights[u]));
        COAM_ADD(scalar_degree_checks,c.length);
        return d;
    }
    unsigned ordered_degree(const int32_t* row) const {
        for (unsigned u : c.order) {
            COAM_ADD(scalar_degree_checks,1);
            if (row[u]) return c.weights[u];
        }
        return 0;
    }

    uint64_t evaluate() {
        make_coordinates();
        if (method == Method::ReferenceProposed) {
            if (metric != Metric::Minimum) throw std::runtime_error("use apply for incremental state");
            std::copy(coord.begin(),coord.end(),work.begin());
            std::array<unsigned,16> rows{};
            std::iota(rows.begin(),rows.end(),0);
            unsigned rank = 0;
            for (unsigned u : c.order) {
                unsigned pivot = rank;
                while (pivot < c.m) {
                    COAM_ADD(pivot_checks,1);
                    if (work[rows[pivot]*c.length+u]) break;
                    ++pivot;
                }
                if (pivot == c.m) continue;
                std::swap(rows[rank],rows[pivot]);
                const auto* src = work.data()+rows[rank]*c.length;
                for (unsigned r = rank+1; r < c.m; ++r) {
                    auto* dst = work.data()+rows[r]*c.length;
                    COAM_ADD(pivot_checks,1);
                    if (dst[u]) {
                        for (unsigned j = 0; j < c.length; ++j) dst[j] ^= src[j];
                        COAM_ADD(scalar_elimination_xor,c.length);
                    }
                }
                if (++rank == c.m) return c.weights[u];
            }
            return 0; // Includes all rank-deficient cases; optional lambda output excluded.
        }
        if (metric == Metric::Maximum) {
            unsigned d = 0;
            for (unsigned b = 0; b < c.m; ++b)
                d = std::max(d,full_degree(coord.data()+b*c.length));
            return d;
        }
        unsigned minimum = c.n;
        hist.fill(0);
        for (unsigned l = 1; l < c.components; ++l) {
            std::fill(combination.begin(),combination.end(),0);
            for (unsigned selected = l; selected; selected &= selected-1) {
                const auto* row = coord.data()+std::countr_zero(selected)*c.length;
                for (unsigned u = 0; u < c.length; ++u) combination[u] ^= row[u];
                COAM_ADD(scalar_component_xor,c.length);
            }
            unsigned d = full_degree(combination.data());
            if (metric == Metric::Spectrum) ++hist[d];
            else minimum = std::min(minimum,d);
        }
        return metric == Metric::Spectrum ? spectrum_sum(hist) : minimum;
    }

    void make_delta(unsigned p, unsigned q) {
        delta.clear();
        top_delta.clear();
        // Enumerate each influence set, excluding its intersection with the other.
        // This avoids a dimension-dependent 3^n precomputed incidence table.
        for (auto [a,b] : {std::pair{p,q},std::pair{q,p}})
            for (unsigned u = a; u < c.length; u = (u+1)|a) {
                COAM_ADD(delta_candidates,1);
                if ((u & b) != b) delta.push_back(u);
            }
        for (unsigned v = p ^ q; v; v &= v-1)
            top_delta.push_back((c.length-1) ^ (1U << std::countr_zero(v)));
        COAM_ADD(delta_coefficients,delta.size());
    }
    void xor_delta(int32_t* row, bool rollback = false) const {
        for (unsigned u : delta) row[u] ^= 1;
#ifdef COAM_COUNT
        if (rollback) COAM_ADD(scalar_rollback_xor,delta.size());
        else COAM_ADD(scalar_coefficient_updates,delta.size());
#else
        (void)rollback;
#endif
    }
    unsigned update_row(int32_t* row, unsigned before) const {
        xor_delta(row);
        const unsigned affected = c.n-1;
        if (before < affected) return affected;
        if (before > affected) return before;
        for (unsigned u : top_delta) {
            COAM_ADD(scalar_guard_checks,1);
            if (row[u]) { COAM_ADD(scalar_guard_kept,1); return before; }
        }
        COAM_ADD(scalar_guard_fallback,1);
        unsigned after = ordered_degree(row);
        if (after < before) COAM_ADD(scalar_degree_drops,1);
        return after;
    }
    void update_truth(unsigned p, unsigned q) {
        for (unsigned b = 0; b < c.m; ++b) if ((delta_output >> b) & 1U) {
            truth[b*c.length+p] ^= 1;
            truth[b*c.length+q] ^= 1;
        }
    }
    uint64_t apply(unsigned p, unsigned q, bool reversible) {
        old_value = value;
        if (reversible && metric == Metric::Spectrum) old_hist = hist;
        delta_output = box[p] ^ box[q];
        std::swap(box[p],box[q]);
        if (method == Method::Reference || metric == Metric::Minimum) {
            update_truth(p,q);
            return value = evaluate();
        }
        if (reversible) old_degrees = degrees;
        if (!delta_output) return value;
        make_delta(p,q);
        if (metric == Metric::Maximum) {
            value = 0;
            for (unsigned b = 0; b < c.m; ++b) {
                if ((delta_output >> b) & 1U)
                    degrees[b] = update_row(coord.data()+b*c.length,degrees[b]);
                value = std::max(value,uint64_t(degrees[b]));
            }
        } else {
            for (unsigned l = 1; l < c.components; ++l) {
                if (!(std::popcount(l & delta_output) & 1U)) continue;
                unsigned before = degrees[l];
                unsigned after = update_row(component_anfs.data()+size_t(l)*c.length,before);
                if (before != after) { --hist[before]; ++hist[after]; degrees[l] = after; }
            }
            value = spectrum_sum(hist);
        }
        return value;
    }
    void undo(unsigned p, unsigned q) {
        std::swap(box[p],box[q]);
        if (method == Method::Reference || metric == Metric::Minimum) update_truth(p,q);
        else {
            if (delta_output) {
                if (metric == Metric::Maximum) {
                    for (unsigned b = 0; b < c.m; ++b) if ((delta_output >> b) & 1U)
                        xor_delta(coord.data()+b*c.length,true);
                } else {
                    for (unsigned l = 1; l < c.components; ++l)
                        if (std::popcount(l & delta_output) & 1U)
                            xor_delta(component_anfs.data()+size_t(l)*c.length,true);
                }
            }
            degrees.swap(old_degrees);
        }
        if (metric == Metric::Spectrum) hist = old_hist;
        value = old_value;
    }
};
}
