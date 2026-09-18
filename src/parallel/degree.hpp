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
#include "parallel.hpp"

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
        parallel_for(0, m, uint64_t(m)*length, [&](unsigned b) {
            Word* row = a.data() + b * words;
            for (unsigned x = 0; x < length; ++x)
                row[x / 64] |= Word((s[x] >> b) & 1U) << (x % 64);
        });
    }

    void coordinates(const Box& s, std::vector<Word>& a) const {
        truth_tables(s,a);
        parallel_for(0, m, uint64_t(m)*words*n, [&](unsigned b) { mobius(a.data() + b*words); });
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
    std::vector<std::vector<Word>> scratch;
    std::vector<Partial> partials;
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
        partials.resize(omp_get_max_threads());
        if (method == Method::Bitwise && metric != Metric::Maximum)
            scratch.assign(omp_get_max_threads(), std::vector<Word>(c.words));
        delta.reserve(c.words);
        if (packed_input) { truth.resize(c.m*c.words); c.truth_tables(box,truth); }
        if (method == Method::Proposed && metric != Metric::Minimum) {
            c.coordinates(box, coord);
            if (metric == Metric::Maximum) {
                degrees.resize(c.m);
                for (unsigned b = 0; b < c.m; ++b) {
                    degrees[b] = uint8_t(c.degree(coord.data() + b * c.words));
                    value = std::max(value, uint64_t(degrees[b]));
                }
            } else {
                // Store only component ANFs. Truth values are parity(lambda & box[x]).
                component_anfs.resize(size_t(c.components) * c.words);
                degrees.resize(c.components);
                // A dependency barrier separates each doubling stage. The
                // highest-bit parent always belongs to a completed stage.
                for (unsigned base = 1, bit = 0; base < c.components; base *= 2, ++bit) {
                    parallel_for(base, 2*base, uint64_t(base)*c.words, [&](unsigned l) {
                        Word* dst = component_anfs.data() + size_t(l)*c.words;
                        const Word* src = component_anfs.data() + size_t(l-base)*c.words;
                        const Word* row = coord.data() + bit*c.words;
                        for (size_t w = 0; w < c.words; ++w) dst[w] = src[w] ^ row[w];
                        degrees[l] = uint8_t(c.degree(dst));
                    });
                }
                for (unsigned l = 1; l < c.components; ++l) ++hist[degrees[l]];
                value = spectrum_sum(hist);
            }
            old_degrees.resize(degrees.size());
        } else value = evaluate();
    }

    void make_coordinates() {
        if (!packed_input) { c.coordinates(box,coord); return; }
        std::copy(truth.begin(),truth.end(),coord.begin());
        parallel_for(0, c.m, uint64_t(c.m)*c.words*c.n, [&](unsigned b) { c.mobius(coord.data()+b*c.words); });
    }

    void update_truth(unsigned p, unsigned q) {
        if (!packed_input) return;
        for (unsigned b = 0; b < c.m; ++b) if ((delta_output >> b) & 1U) {
            truth[b*c.words+p/64] ^= Word(1) << (p%64);
            truth[b*c.words+q/64] ^= Word(1) << (q%64);
        }
    }

    // Algorithm 2, numerical mode: reverse Gaussian elimination of packed rows.
    // The optional lambda enumeration below is deliberately outside this path.
    unsigned minimum() {
        make_coordinates();
        std::copy(coord.begin(), coord.end(), work.begin());
        std::array<unsigned, 16> rows{};
        std::iota(rows.begin(), rows.end(), 0);
        unsigned rank = 0;
        for (unsigned u : c.order) {
            unsigned pivot = rank;
            while (pivot < c.m && !((work[rows[pivot]*c.words + u/64] >> (u%64)) & 1U)) ++pivot;
            if (pivot == c.m) continue;
            std::swap(rows[rank], rows[pivot]);
            const Word* src = work.data() + rows[rank] * c.words;
            parallel_for(rank+1, c.m, uint64_t(c.m-rank-1)*c.words, [&](unsigned r) {
                Word* dst = work.data() + rows[r] * c.words;
                if ((dst[u/64] >> (u%64)) & 1U) {
                    for (size_t w = 0; w < c.words; ++w) dst[w] ^= src[w];
                    COAM_ADD(word_elimination_xor,c.words);
                }
            });
            if (++rank == c.m) return c.weights[u];
        }
        // A nonzero lambda with zero ANF exists whenever rank < m.
        return 0;
    }

    uint64_t bitwise() {
        make_coordinates();
        if (metric == Metric::Maximum) {
            unsigned d = 0;
            for (unsigned b = 0; b < c.m; ++b)
                d = std::max(d, c.degree(coord.data() + b*c.words));
            return d;
        }
        hist.fill(0);
        unsigned blocks = parallel_blocks(1, c.components, uint64_t(c.components-1)*c.words,
            [&](unsigned first, unsigned last, unsigned block) {
                auto& buffer = scratch[block];
                auto& part = partials[block];
                std::fill(buffer.begin(), buffer.end(), 0);
                part.histogram.fill(0);
                part.minimum = c.n;
                // Seed the predecessor of this Gray segment, then retain
                // one coordinate XOR per successive component within it.
                unsigned predecessor = first-1;
                unsigned gray = predecessor ^ (predecessor >> 1);
                for (unsigned mask = gray; mask; mask &= mask-1) {
                    const Word* row = coord.data() + std::countr_zero(mask)*c.words;
                    for (size_t w = 0; w < c.words; ++w) buffer[w] ^= row[w];
                }
                for (unsigned i = first; i < last; ++i) {
                    const Word* row = coord.data() + std::countr_zero(i)*c.words;
                    for (size_t w = 0; w < c.words; ++w) buffer[w] ^= row[w];
                    unsigned d = c.degree(buffer.data());
                    if (metric == Metric::Spectrum) ++part.histogram[d];
                    else {
                        part.minimum = std::min(part.minimum, d);
                        if (!part.minimum) break;
                    }
                }
            });
        unsigned min_degree = c.n;
        for (unsigned block = 0; block < blocks; ++block) {
            min_degree = std::min(min_degree, partials[block].minimum);
            for (unsigned d = 0; d <= c.n; ++d) hist[d] += partials[block].histogram[d];
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
        for (size_t w = 0; w < c.words; ++w) {
            Word a = ((w & (p >> 6)) == (p >> 6)) ? c.supermask[p & 63] : 0;
            Word b = ((w & (q >> 6)) == (q >> 6)) ? c.supermask[q & 63] : 0;
            Word x = a ^ b;
            if (c.n < 6) x &= (Word(1) << c.length) - 1;
            if (x) delta.emplace_back(w, x);
            COAM_ADD(delta_coefficients,std::popcount(x));
            COAM_ADD(delta_candidates,1);
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
        if (reversible) old_degrees = degrees;
        if (!delta_output) return value;
        make_delta(p, q);
        if (metric == Metric::Maximum) {
            value = 0;
            for (unsigned b = 0; b < c.m; ++b) {
                if ((delta_output >> b) & 1U) {
                    Word* row = coord.data() + b * c.words;
                    xor_delta(row);
                    degrees[b] = uint8_t(updated_degree(row,degrees[b]));
                }
                value = std::max(value, uint64_t(degrees[b]));
            }
        } else {
            unsigned blocks = parallel_blocks(1, c.components,
                uint64_t(c.components-1)*std::max(size_t(1), delta.size()),
                [&](unsigned first, unsigned last, unsigned block) {
                    auto& changes = partials[block].histogram;
                    changes.fill(0);
                    for (unsigned l = first; l < last; ++l) {
                        if (!(std::popcount(l & delta_output) & 1U)) continue;
                        Word* row = component_anfs.data() + size_t(l)*c.words;
                        xor_delta(row);
                        unsigned before = degrees[l], after = updated_degree(row,before);
                        if (before != after) {
                            --changes[before]; ++changes[after]; degrees[l] = uint8_t(after);
                        }
                    }
                });
            for (unsigned block = 0; block < blocks; ++block)
                for (unsigned d = 0; d <= c.n; ++d)
                    hist[d] = uint32_t(int64_t(hist[d]) + partials[block].histogram[d]);
            value = spectrum_sum(hist);
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
                    parallel_for(1, c.components, uint64_t(c.components-1)*delta.size(), [&](unsigned l) {
                        if (std::popcount(l & delta_output) & 1U)
                            xor_delta(component_anfs.data() + size_t(l)*c.words,true);
                    });
                }
            }
            degrees.swap(old_degrees);
        }
        if (metric == Metric::Spectrum) hist = old_hist;
        value = old_value;
    }

    // Algorithm 2's optional output: enumerate ker(A_{>minimum}) \ {0}.
    std::vector<unsigned> minimum_lambdas() {
        unsigned d = minimum();
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
