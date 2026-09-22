#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "counts.hpp"

namespace coam {
using Word = uint64_t;
using Box = std::vector<uint32_t>;
using Hist = std::array<uint32_t, 17>;
enum class Metric { Minimum, Maximum, Spectrum };
enum class Method { Proposed, Bitwise, Peigen, Reference, ReferenceProposed };

inline Metric metric_of(const std::string& s) {
    if (s == "min") return Metric::Minimum;
    if (s == "max") return Metric::Maximum;
    if (s == "spectrum") return Metric::Spectrum;
    throw std::runtime_error("unknown metric: " + s);
}
inline Method method_of(const std::string& s) {
    if (s == "proposed") return Method::Proposed;
    if (s == "bitwise") return Method::Bitwise;
    if (s == "peigen") return Method::Peigen;
    if (s == "reference") return Method::Reference;
    if (s == "reference-proposed") return Method::ReferenceProposed;
    throw std::runtime_error("unknown method: " + s);
}
inline uint64_t spectrum_sum(const Hist& h) {
    uint64_t s = 0;
    for (unsigned d = 0; d < h.size(); ++d) s += uint64_t(d) * h[d];
    return s;
}
inline unsigned hist_min(const Hist& h) {
    for (unsigned d = 0; d < h.size(); ++d) if (h[d]) return d;
    throw std::runtime_error("empty spectrum");
}
inline unsigned hist_max(const Hist& h) {
    for (int d = 16; d >= 0; --d) if (h[d]) return unsigned(d);
    throw std::runtime_error("empty spectrum");
}

// Dimension-only data. Construction is recorded separately from evaluation.
struct Context {
    unsigned n, m, length, components;
    size_t words;
    std::vector<unsigned> order;
    std::vector<uint8_t> weights;
    std::vector<Word> degree_masks;
    std::array<Word, 64> supermask{};
    Context(unsigned ni, unsigned mi) : n(ni), m(mi) {
        if (n < 1 || n > 16 || m < 1 || m > 16)
            throw std::runtime_error("supported dimensions: 1..16");
        length = 1U << n;
        components = 1U << m;
        words = (length + 63) / 64;
        weights.resize(length);
        degree_masks.assign((n + 1) * words, 0);
        std::array<std::vector<unsigned>, 17> layers;
        for (unsigned u = 0; u < length; ++u) {
            unsigned d = std::popcount(u);
            weights[u] = uint8_t(d);
            layers[d].push_back(u);
            degree_masks[d * words + u / 64] |= Word(1) << (u % 64);
        }
        order.reserve(length);
        for (int d = int(n); d >= 0; --d)
            order.insert(order.end(), layers[d].rbegin(), layers[d].rend());
        for (unsigned p = 0; p < 64; ++p)
            for (unsigned u = 0; u < 64; ++u)
                if ((u & p) == p) supermask[p] |= Word(1) << u;
    }

    // Packed Boolean Mobius transform: six within-word stages, then word XORs.
    void mobius(Word* a) const {
        COAM_ADD(word_mobius_xor,words*std::min(n,6U) + (n > 6 ? (n-6)*(words/2) : 0));
        constexpr Word masks[] = {0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL,
            0xf0f0f0f0f0f0f0f0ULL, 0xff00ff00ff00ff00ULL,
            0xffff0000ffff0000ULL, 0xffffffff00000000ULL};
        for (size_t w = 0; w < words; ++w) {
            Word x = a[w];
            for (unsigned b = 0; b < std::min(n, 6U); ++b)
                x ^= (x << (1U << b)) & masks[b];
            a[w] = x;
        }
        for (size_t step = 1; step < words; step *= 2)
            for (size_t base = 0; base < words; base += 2 * step)
                for (size_t j = 0; j < step; ++j)
                    a[base + step + j] ^= a[base + j];
    }

    void truth_tables(const Box& s, std::vector<Word>& a) const {
        std::fill(a.begin(), a.end(), 0);
        for (unsigned b = 0; b < m; ++b) {
            Word* row = a.data() + b * words;
            for (unsigned x = 0; x < length; ++x)
                row[x / 64] |= Word((s[x] >> b) & 1U) << (x % 64);
        }
    }

    void coordinates(const Box& s, std::vector<Word>& a) const {
        truth_tables(s,a);
        for (unsigned b = 0; b < m; ++b) mobius(a.data() + b*words);
    }

    // Bakoev's bitwise WLO for small vectors; CB-WLO for larger vectors.
    // deg(0)=deg(1)=0, as explicitly adopted by the manuscript.
    unsigned degree(const Word* a) const {
        if (n <= 8) {
            for (unsigned d = n; d > 0; --d)
                for (size_t w = 0; w < words; ++w)
                    { COAM_ADD(word_degree_checks,1); if (a[w] & degree_masks[d * words + w]) return d; }
        } else {
            for (unsigned u : order) {
                COAM_ADD(word_degree_checks,1);
                if ((a[u / 64] >> (u % 64)) & 1U) return weights[u];
            }
        }
        return 0;
    }
};

// Defined in a separate translation unit containing unmodified upstream headers.
uint64_t peigen_evaluate(const Box&, unsigned n, Metric, Hist&);
uint64_t peigen_prepared(const std::vector<Word>&, unsigned n, Metric, Hist&);
Hist peigen_native_spectrum(const Box&, unsigned n);

struct Engine {
    const Context& c;
    Method method;
    Metric metric;
    bool packed_input;
    Box box;
    std::vector<Word> coord, work, component_anfs, combination, truth;
    std::vector<uint8_t> degrees, old_degrees;
    std::vector<std::pair<unsigned, uint8_t>> changed_degrees;
    // Four 16-bit output coefficient vectors per word (Algorithm 2).
    std::vector<Word> output_anfs, output_truth;
    // A tile stores 64 labelled components, one word per ANF monomial.
    std::vector<Word> degree_lanes;
    std::vector<unsigned> delta_terms;
    unsigned component_groups = 0;
    std::vector<std::pair<size_t, Word>> delta;
    std::array<std::pair<size_t, Word>,16> top_delta{};
    unsigned top_count = 0;
    Hist hist{}, old_hist{};
    uint64_t value = 0, old_value = 0;
    unsigned delta_output = 0;

    Engine(const Context& cx, Method me, Metric mt, const Box& s, bool prepared = false)
        : c(cx), method(me), metric(mt),
          packed_input(prepared && !(me == Method::Proposed && mt != Metric::Minimum)),
          box(s), coord(c.m * c.words),
          work(c.m * c.words), combination(c.words) {
        if (box.size() != c.length ||
            std::any_of(box.begin(), box.end(), [&](uint32_t x) {return x >= c.components;}))
            throw std::runtime_error("invalid truth table");
        if (method == Method::Peigen && (c.n < 3 || c.n > 8 || c.n != c.m))
            throw std::runtime_error("PEIGEN requires square 3..8-bit input");
        delta.reserve(c.words);
        if (packed_input) { truth.resize(c.m*c.words); c.truth_tables(box,truth); }
#ifdef COAM_OUTPUT_PACKED_MINIMUM
        if (method == Method::Proposed && metric == Metric::Minimum) {
            output_anfs.resize((c.length + 3) / 4);
            if (packed_input) {
                output_truth.resize(output_anfs.size());
                pack_outputs(output_truth);
            }
        }
#endif
        if (method == Method::Proposed && metric != Metric::Minimum) {
            c.coordinates(box, coord);
            if (metric == Metric::Maximum) {
                degrees.resize(c.m);
                for (unsigned b = 0; b < c.m; ++b) {
                    degrees[b] = uint8_t(c.degree(coord.data() + b * c.words));
                    value = std::max(value, uint64_t(degrees[b]));
                }
            } else {
#ifdef COAM_TILED_SPECTRUM
                if (c.m >= 6) initialize_tiles();
                else
#endif
                {
                // Store only component ANFs. Truth values are parity(lambda & box[x]).
                component_anfs.resize(size_t(c.components) * c.words);
                degrees.resize(c.components);
                for (unsigned l = 1; l < c.components; ++l) {
                    unsigned bit = std::countr_zero(l), parent = l & (l - 1);
                    Word* dst = component_anfs.data() + size_t(l) * c.words;
                    const Word* src = component_anfs.data() + size_t(parent) * c.words;
                    const Word* row = coord.data() + bit * c.words;
                    for (size_t w = 0; w < c.words; ++w) dst[w] = src[w] ^ row[w];
                    COAM_ADD(word_component_xor,c.words);
                    degrees[l] = uint8_t(c.degree(dst));
                    ++hist[degrees[l]];
                }
                value = spectrum_sum(hist);
                }
            }
            old_degrees.resize(degrees.size());
            if (metric == Metric::Spectrum) changed_degrees.reserve(c.components / 2);
        } else value = evaluate();
    }

    void make_coordinates() {
#ifdef COAM_OUTPUT_PACKED_MINIMUM
        if (method == Method::Proposed && metric == Metric::Minimum) {
            c.coordinates(box, coord);
            return;
        }
#endif
        if (!packed_input) { c.coordinates(box,coord); return; }
        std::copy(truth.begin(),truth.end(),coord.begin());
        for (unsigned b = 0; b < c.m; ++b) c.mobius(coord.data()+b*c.words);
    }

    void update_truth(unsigned p, unsigned q) {
        if (!packed_input) return;
#ifdef COAM_OUTPUT_PACKED_MINIMUM
        if (method == Method::Proposed && metric == Metric::Minimum) {
            output_truth[p / 4] ^= Word(delta_output) << (16 * (p % 4));
            output_truth[q / 4] ^= Word(delta_output) << (16 * (q % 4));
            return;
        }
#endif
        for (unsigned b = 0; b < c.m; ++b) if ((delta_output >> b) & 1U) {
            truth[b*c.words+p/64] ^= Word(1) << (p%64);
            truth[b*c.words+q/64] ^= Word(1) << (q%64);
        }
    }

    // Algorithm 2, numerical mode: reverse Gaussian elimination of packed rows.
    // The optional lambda enumeration below is deliberately outside this path.
    unsigned minimum(bool preserve_coordinates = false) {
#ifdef COAM_OUTPUT_PACKED_MINIMUM
        if (!preserve_coordinates) return output_packed_minimum();
#endif
        make_coordinates();
        if (preserve_coordinates) std::copy(coord.begin(), coord.end(), work.begin());
        else coord.swap(work);
        std::array<unsigned, 16> rows{};
        std::iota(rows.begin(), rows.end(), 0);
        unsigned rank = 0;
        for (unsigned u : c.order) {
            unsigned pivot = rank;
            while (pivot < c.m && !((work[rows[pivot]*c.words + u/64] >> (u%64)) & 1U)) ++pivot;
            if (pivot == c.m) continue;
            std::swap(rows[rank], rows[pivot]);
            const Word* src = work.data() + rows[rank] * c.words;
            for (unsigned r = rank + 1; r < c.m; ++r) {
                Word* dst = work.data() + rows[r] * c.words;
                if ((dst[u/64] >> (u%64)) & 1U) {
                    for (size_t w = 0; w < c.words; ++w) dst[w] ^= src[w];
                    COAM_ADD(word_elimination_xor,c.words);
                }
            }
            if (++rank == c.m) return c.weights[u];
        }
        // A nonzero lambda with zero ANF exists whenever rank < m.
        return 0;
    }

    void pack_outputs(std::vector<Word>& dst) const {
        if (c.length < 4) {
            dst[0] = Word(box[0]) | (Word(box[1]) << 16);
            return;
        }
        for (unsigned u = 0; u < c.length; u += 4)
            dst[u / 4] = Word(box[u]) | (Word(box[u+1]) << 16) |
                         (Word(box[u+2]) << 32) | (Word(box[u+3]) << 48);
    }

    unsigned output_packed_minimum() {
        if (packed_input) std::copy(output_truth.begin(), output_truth.end(), output_anfs.begin());
        else pack_outputs(output_anfs);
        const size_t size = output_anfs.size();
        for (Word& v : output_anfs) {
            v ^= (v << 16) & 0xffff0000ffff0000ULL;
            if (c.n >= 2) v ^= (v << 32) & 0xffffffff00000000ULL;
        }
        for (size_t step = 1; step < size; step *= 2)
            for (size_t base = 0; base < size; base += 2 * step)
                for (size_t j = 0; j < step; ++j)
                    output_anfs[base + step + j] ^= output_anfs[base + j];
        // Reverse Gaussian elimination, with coefficient columns represented
        // by m-bit integers. The last independent column fixes the minimum.
        std::array<unsigned,16> basis{};
        unsigned rank = 0;
        for (unsigned u : c.order) {
            unsigned v = unsigned(output_anfs[u/4] >> (16*(u%4))) & 65535U;
            while (v) {
                unsigned pivot = 31U - std::countl_zero(v);
                if (basis[pivot]) v ^= basis[pivot];
                else {
                    basis[pivot] = v;
                    if (++rank == c.m) return c.weights[u];
                    break;
                }
            }
        }
        return 0;
    }

    static Word component_mask(unsigned v, unsigned group) {
        static constexpr auto table = [] {
            constexpr std::array<Word,6> low = {0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL,
                0xf0f0f0f0f0f0f0f0ULL, 0xff00ff00ff00ff00ULL,
                0xffff0000ffff0000ULL, 0xffffffff00000000ULL};
            std::array<Word,64> a{};
            for (unsigned i = 1; i < 64; ++i) a[i] = a[i & (i-1)] ^ low[std::countr_zero(i)];
            return a;
        }();
        return table[v & 63U] ^ (Word(0) - Word(std::popcount((v >> 6) & group) & 1U));
    }

    Word valid_lanes(unsigned group) const {
        Word lanes = c.m < 6 ? (Word(1) << c.components) - 1 : ~Word(0);
        return group == 0 ? lanes & ~Word(1) : lanes;
    }

    void initialize_tiles() {
        component_groups = (c.components + 63) / 64;
        component_anfs.resize(size_t(component_groups) * c.length);
        degrees.assign(c.components, 0);
        degree_lanes.assign(size_t(c.n + 1) * component_groups, 0);
        delta_terms.reserve(c.length);
        std::vector<unsigned> coefficients(c.length);
        for (unsigned u = 0; u < c.length; ++u)
            for (unsigned b = 0; b < c.m; ++b)
                coefficients[u] |= unsigned((coord[b*c.words + u/64] >> (u%64)) & 1U) << b;
        for (unsigned g = 0; g < component_groups; ++g) {
            Word* tile = component_anfs.data() + size_t(g) * c.length;
            for (unsigned u = 0; u < c.length; ++u) tile[u] = component_mask(coefficients[u], g);
            Word remaining = valid_lanes(g);
            for (unsigned u : c.order) {
                Word found = tile[u] & remaining;
                unsigned d = c.weights[u];
                degree_lanes[d*component_groups+g] |= found;
                hist[d] += std::popcount(found);
                for (Word bits = found; bits; bits &= bits-1)
                    degrees[g*64 + std::countr_zero(bits)] = uint8_t(d);
                remaining &= ~found;
                if (!remaining) break;
            }
            // Zero components share the manuscript convention deg(0)=0.
            degree_lanes[g] |= remaining;
            hist[0] += std::popcount(remaining);
        }
        value = spectrum_sum(hist);
    }

    void set_lane_degrees(unsigned group, Word bits, unsigned after, bool reversible) {
        while (bits) {
            unsigned bit = std::countr_zero(bits), l = group*64 + bit;
            bits &= bits-1;
            unsigned before = degrees[l];
            if (before == after) continue;
            if (reversible) changed_degrees.emplace_back(l, uint8_t(before));
            Word flag = Word(1) << bit;
            degree_lanes[before*component_groups+group] &= ~flag;
            degree_lanes[after*component_groups+group] |= flag;
            degrees[l] = uint8_t(after);
            --hist[before]; ++hist[after];
            value = value - before + after;
        }
    }

    void apply_tiles(bool reversible) {
        delta_terms.clear();
        for (auto [w, bits] : delta)
            while (bits) {
                delta_terms.push_back(unsigned(w*64 + std::countr_zero(bits)));
                bits &= bits-1;
            }
        for (unsigned g = 0; g < component_groups; ++g) {
            Word mask = component_mask(delta_output, g) & valid_lanes(g);
            if (!mask) continue;
            Word* tile = component_anfs.data() + size_t(g)*c.length;
            for (unsigned u : delta_terms) tile[u] ^= mask;
            const unsigned top = c.n - 1;
            Word at_top = mask & degree_lanes[top*component_groups+g];
            Word below = mask & ~(degree_lanes[c.n*component_groups+g] | at_top);
            set_lane_degrees(g, below, top, reversible);
            Word guard = 0;
            for (unsigned v = pending_input_difference; v; v &= v-1)
                guard |= tile[(c.length-1) ^ (1U << std::countr_zero(v))];
            Word remaining = at_top & ~guard;
            if (!remaining) continue;
            for (unsigned u : c.order) {
                if (c.weights[u] > top) continue;
                Word found = tile[u] & remaining;
                set_lane_degrees(g, found, c.weights[u], reversible);
                remaining &= ~found;
                if (!remaining) break;
            }
            set_lane_degrees(g, remaining, 0, reversible);
        }
    }

    void apply_compact(bool reversible) {
        Word affected = component_mask(delta_output, 0) & valid_lanes(0);
        Word change = delta[0].second, guard = top_delta[0].second;
        while (affected) {
            unsigned l = std::countr_zero(affected);
            affected &= affected-1;
            Word& row = component_anfs[l];
            row ^= change;
            unsigned before = degrees[l], after = before;
            if (before < c.n-1) after = c.n-1;
            else if (before == c.n-1 && !(row & guard)) after = c.degree(&row);
            if (before != after) {
                if (reversible) changed_degrees.emplace_back(l, uint8_t(before));
                degrees[l] = uint8_t(after);
                --hist[before]; ++hist[after];
                value = value - before + after;
            }
        }
    }

    unsigned pending_input_difference = 0;

    bool component_coefficient(unsigned l, unsigned u) const {
#ifdef COAM_TILED_SPECTRUM
        if (component_groups)
            return (component_anfs[size_t(l/64)*c.length+u] >> (l%64)) & 1U;
#endif
        return (component_anfs[size_t(l)*c.words+u/64] >> (u%64)) & 1U;
    }

    uint64_t bitwise() {
        make_coordinates();
        if (metric == Metric::Maximum) {
            unsigned d = 0;
            for (unsigned b = 0; b < c.m; ++b)
                d = std::max(d, c.degree(coord.data() + b*c.words));
            return d;
        }
        std::fill(combination.begin(), combination.end(), 0);
        hist.fill(0);
        unsigned min_degree = c.n;
        // Gray-code enumeration: one packed coordinate XOR per nonzero component.
        for (unsigned i = 1; i < c.components; ++i) {
            const Word* row = coord.data() + std::countr_zero(i) * c.words;
            for (size_t w = 0; w < c.words; ++w) combination[w] ^= row[w];
            COAM_ADD(word_component_xor,c.words);
            unsigned d = c.degree(combination.data());
            if (metric == Metric::Spectrum) ++hist[d];
            else {
                min_degree = std::min(min_degree, d);
                if (min_degree == 0) break;
            }
        }
        return metric == Metric::Spectrum ? spectrum_sum(hist) : min_degree;
    }



    uint64_t evaluate() {
        if (method == Method::Peigen)
            return packed_input ? peigen_prepared(truth,c.n,metric,hist) : peigen_evaluate(box,c.n,metric,hist);
        if (method == Method::Bitwise) return bitwise();
        if (metric == Metric::Minimum) return minimum();
        throw std::runtime_error("incremental state must be updated with apply");
    }

    // A(p) triangle A(q), represented as a sparse list of changed 64-bit words.
    void make_delta(unsigned p, unsigned q) {
        delta.clear();
        top_count = 0;
        pending_input_difference = p ^ q;
        auto append = [&](size_t w) {
            Word a = ((w & (p >> 6)) == (p >> 6)) ? c.supermask[p & 63] : 0;
            Word b = ((w & (q >> 6)) == (q >> 6)) ? c.supermask[q & 63] : 0;
            Word x = a ^ b;
            if (c.n < 6) x &= (Word(1) << c.length) - 1;
            if (x) delta.emplace_back(w, x);
            COAM_ADD(delta_coefficients,std::popcount(x));
            COAM_ADD(delta_candidates,1);
        };
#ifdef COAM_SPARSE_DELTA
        unsigned ph = p >> 6, qh = q >> 6;
        if (c.words >= 16 && (c.words >> std::popcount(ph)) +
                              (c.words >> std::popcount(qh)) < c.words) {
            for (unsigned w = ph; w < c.words; w = (w+1) | ph) append(w);
            for (unsigned w = qh; w < c.words; w = (w+1) | qh)
                if ((w & ph) != ph) append(w);
        } else
#endif
        for (size_t w = 0; w < c.words; ++w) {
            append(w);
        }
        // Weight n-1 members of A(p) triangle A(q) omit one bit on which p,q differ.
        for (unsigned v = p ^ q; v; v &= v-1) {
            unsigned u = (c.length-1) ^ (1U << std::countr_zero(v));
            size_t w = u/64;
            Word mask = Word(1) << (u%64);
            if (top_count && top_delta[top_count-1].first == w) top_delta[top_count-1].second |= mask;
            else top_delta[top_count++] = {w,mask};
        }
    }

    void xor_delta(Word* row, bool rollback = false) const {
        for (auto [w, bits] : delta) row[w] ^= bits;
#ifdef COAM_COUNT
        if (rollback) COAM_ADD(word_rollback_xor,delta.size());
        else COAM_ADD(word_coefficient_updates,delta.size());
#else
        (void)rollback;
#endif
    }

    bool highest_affected_survives(const Word* row) const {
        for (unsigned i = 0; i < top_count; ++i) {
            auto [w,mask] = top_delta[i];
            COAM_ADD(top_guard_checks,1);
            if (row[w] & mask) { COAM_ADD(top_guard_kept,1); return true; }
        }
        COAM_ADD(top_guard_fallback,1);
        return false;
    }

    unsigned updated_degree(const Word* row, unsigned before) const {
        const unsigned affected = c.n-1;
        if (before < affected) return affected;
        if (before > affected || highest_affected_survives(row)) return before;
        unsigned after = c.degree(row);
        if (after < before) COAM_ADD(degree_drops,1);
        return after;
    }

    uint64_t apply(unsigned p, unsigned q, bool reversible) {
        old_value = value;
        if (reversible && metric == Metric::Spectrum) old_hist = hist;
        delta_output = box[p] ^ box[q];
        std::swap(box[p], box[q]);
        update_truth(p,q);
        if (method != Method::Proposed || metric == Metric::Minimum)
            return value = evaluate();
        if (reversible) {
            if (metric == Metric::Maximum) old_degrees = degrees;
            else changed_degrees.clear();
        }
        if (!delta_output) return value;
        if (metric == Metric::Maximum && c.words == 1 && c.m <= 4) {
            // Tiny coordinate sets fit in one word: fuse update and maximum.
            Word change = c.supermask[p] ^ c.supermask[q];
            if (c.n < 6) change &= (Word(1) << c.length) - 1;
            delta.clear();
            delta.emplace_back(0, change);
            const unsigned top = c.n - 1;
            value = 0;
            for (unsigned b = 0; b < c.m; ++b) {
                if ((delta_output >> b) & 1U) {
                    coord[b] ^= change;
                    unsigned before = degrees[b];
                    if (before < top) degrees[b] = uint8_t(top);
                    else if (before == top && !(coord[b] & c.degree_masks[top]))
                        degrees[b] = uint8_t(c.degree(coord.data() + b));
                }
                value = std::max(value, uint64_t(degrees[b]));
            }
            return value;
        }
        make_delta(p, q);
        if (metric == Metric::Maximum) {
            bool dropped = false;
            for (unsigned bits = delta_output; bits; bits &= bits-1) {
                unsigned b = std::countr_zero(bits), before = degrees[b];
                Word* row = coord.data() + b * c.words;
                xor_delta(row);
                unsigned after = updated_degree(row, before);
                degrees[b] = uint8_t(after);
                value = std::max(value, uint64_t(after));
                dropped |= before == old_value && after < before;
            }
            if (dropped && value == old_value)
                value = *std::max_element(degrees.begin(), degrees.end());
        } else {
#ifdef COAM_TILED_SPECTRUM
            if (component_groups) apply_tiles(reversible);
            else if (c.words == 1 && c.m < 6) apply_compact(reversible);
            else
#endif
            {
#ifdef COAM_DIRECT_COMPONENTS
            // Solve lambda dot delta_output = 1 using its lowest nonzero bit.
            const unsigned pivot = delta_output & -delta_output;
            for (unsigned j = 0; j < c.components / 2; ++j) {
                unsigned l = (j & (pivot - 1)) | ((j & ~(pivot - 1)) << 1);
                l |= ((std::popcount(l & delta_output) & 1U) ^ 1U) * pivot;
#else
            for (unsigned l = 1; l < c.components; ++l) {
                if ((std::popcount(l & delta_output) & 1U) == 0) continue;
#endif
                Word* row = component_anfs.data() + size_t(l) * c.words;
                xor_delta(row);
                unsigned before = degrees[l], after = updated_degree(row,before);
                if (before != after) {
                    if (reversible) changed_degrees.emplace_back(l, uint8_t(before));
                    --hist[before]; ++hist[after]; degrees[l] = uint8_t(after);
                }
            }
            value = spectrum_sum(hist);
            }
        }
        return value;
    }

    void undo(unsigned p, unsigned q) {
        std::swap(box[p], box[q]);
        update_truth(p,q);
        if (method == Method::Proposed && metric != Metric::Minimum) {
            if (delta_output) {
                if (metric == Metric::Maximum) {
                    for (unsigned b = 0; b < c.m; ++b)
                        if ((delta_output >> b) & 1U) xor_delta(coord.data() + b*c.words,true);
                } else {
#ifdef COAM_TILED_SPECTRUM
                    if (component_groups) {
                    for (unsigned g = 0; g < component_groups; ++g) {
                        Word mask = component_mask(delta_output, g) & valid_lanes(g);
                        if (!mask) continue;
                        Word* tile = component_anfs.data() + size_t(g)*c.length;
                        for (unsigned u : delta_terms) tile[u] ^= mask;
                    }
                    } else if (c.words == 1 && c.m < 6) {
                        Word affected = component_mask(delta_output, 0) & valid_lanes(0);
                        while (affected) {
                            unsigned l = std::countr_zero(affected);
                            affected &= affected-1;
                            component_anfs[l] ^= delta[0].second;
                        }
                    } else
#endif
                    {
#ifdef COAM_DIRECT_COMPONENTS
                    const unsigned pivot = delta_output & -delta_output;
                    for (unsigned j = 0; j < c.components / 2; ++j) {
                        unsigned l = (j & (pivot - 1)) | ((j & ~(pivot - 1)) << 1);
                        l |= ((std::popcount(l & delta_output) & 1U) ^ 1U) * pivot;
                        xor_delta(component_anfs.data() + size_t(l)*c.words,true);
                    }
#else
                    for (unsigned l = 1; l < c.components; ++l)
                        if (std::popcount(l & delta_output) & 1U)
                            xor_delta(component_anfs.data() + size_t(l)*c.words,true);
#endif
                    }
                }
            }
            if (metric == Metric::Maximum) degrees.swap(old_degrees);
            else for (auto [l, before] : changed_degrees) {
#ifdef COAM_TILED_SPECTRUM
                if (component_groups) {
                    Word bit = Word(1) << (l%64);
                    degree_lanes[degrees[l]*component_groups+l/64] &= ~bit;
                    degree_lanes[before*component_groups+l/64] |= bit;
                }
#endif
                degrees[l] = before;
            }
        }
        if (metric == Metric::Spectrum) hist = old_hist;
        value = old_value;
    }

    // Algorithm 2's optional output: enumerate ker(A_{>minimum}) \ {0}.
    std::vector<unsigned> minimum_lambdas() {
        unsigned d = minimum(true);
        std::array<unsigned, 16> basis{};
        for (unsigned u : c.order) {
            if (c.weights[u] <= d) break;
            unsigned v = 0;
            for (unsigned b = 0; b < c.m; ++b)
                v |= unsigned((coord[b*c.words + u/64] >> (u%64)) & 1U) << b;
            while (v) {
                unsigned b = 31U - std::countl_zero(v);
                if (basis[b]) v ^= basis[b];
                else { basis[b] = v; break; }
            }
        }
        std::vector<unsigned> kernel{0};
        for (unsigned f = 0; f < c.m; ++f) if (!basis[f]) {
            unsigned x = 1U << f;
            for (unsigned b = 0; b < c.m; ++b)
                if (basis[b] && (std::popcount(basis[b] & x) & 1U)) x ^= 1U << b;
            size_t old_size = kernel.size();
            for (size_t j = 0; j < old_size; ++j) kernel.push_back(kernel[j] ^ x);
        }
        kernel.erase(kernel.begin());
        std::sort(kernel.begin(), kernel.end());
        return kernel;
    }
};
} // namespace coam
