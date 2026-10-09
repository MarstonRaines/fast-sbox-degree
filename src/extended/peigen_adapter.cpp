#include "degree.hpp"
#include "PEIGEN/func.hpp"

namespace coam {
template<int N>
Peigen::function_t<N> peigen_box(const Box& box) {
    // N=3's upstream constructor reads one whole 16-byte SIMD register.
    alignas(32) std::array<uint8_t, (1 << N) < 16 ? 16 : (1 << N)> lut{};
    std::copy(box.begin(), box.end(), lut.begin());
    return Peigen::function_t<N>(lut.data());
}

template<int N>
Hist native_spectrum(const Box& box) {
    auto s = peigen_box<N>(box);
    int counts[N + 1], maximum, minimum;
    s.degree(counts, maximum, minimum);
    Hist h{};
    for (int d = 0; d <= N; ++d) h[d] = counts[d];
    // Upstream SIMD branches count lambda=0; the other branches exclude it.
    if constexpr (N == 3 || N == 4) --h[0];
    return h;
}

template<int N>
uint64_t evaluate_core(const Peigen::function_t<N>& s, Metric metric, Hist& h) {
    if (metric == Metric::Maximum) {
        bit_slice_t<N> a{};
        s.get_coordinates_ANF(a);
        int d = 0;
        for (auto& row : a) d = std::max(d, s.degree_from_ANF(row));
        return d;
    }
    if constexpr (N <= 4) {
        int counts[N + 1], maximum, minimum;
        s.degree(counts, maximum, minimum);
        if (metric == Metric::Minimum) return minimum;
        h.fill(0);
        for (int d = 0; d <= N; ++d) h[d] = counts[d];
        --h[0];
    } else {
        // Use the library's WLO degree core rather than its slower full scan in
        // degree(). Both calls below are unmodified upstream numerical methods.
        std::array<bit_slice_l_t<N>, 1 << N> a{};
        s.get_components_ANF(a);
        unsigned minimum = N;
        h.fill(0);
        for (unsigned l = 1; l < (1U << N); ++l) {
            unsigned d = s.degree_from_ANF(a[l]);
            if (metric == Metric::Spectrum) ++h[d];
            else minimum = std::min(minimum, d);
        }
        if (metric == Metric::Minimum) return minimum;
    }
    return spectrum_sum(h);
}

template<int N>
uint64_t evaluate(const Box& box, Metric metric, Hist& h) {
    return evaluate_core(peigen_box<N>(box),metric,h);
}

template<int N>
uint64_t prepared(const std::vector<Word>& truth, Metric metric, Hist& h) {
    Peigen::function_t<N> s;
    constexpr size_t words = ((1U << N) + 63) / 64;
    for (unsigned b = 0; b < N; ++b)
        for (size_t w = 0; w < words; ++w)
            s.bit_slice[b][w] = UINT_<N>(truth[b*words+w]);
    return evaluate_core(s,metric,h);
}

uint64_t peigen_evaluate(const Box& box, unsigned n, Metric metric, Hist& h) {
    switch (n) {
#define CASE(N) case N: return evaluate<N>(box, metric, h)
        CASE(3); CASE(4); CASE(5); CASE(6); CASE(7); CASE(8);
#undef CASE
        default: throw std::runtime_error("PEIGEN dimension out of range");
    }
}
uint64_t peigen_prepared(const std::vector<Word>& truth, unsigned n, Metric metric, Hist& h) {
    switch (n) {
#define CASE(N) case N: return prepared<N>(truth, metric, h)
        CASE(3); CASE(4); CASE(5); CASE(6); CASE(7); CASE(8);
#undef CASE
        default: throw std::runtime_error("PEIGEN dimension out of range");
    }
}
Hist peigen_native_spectrum(const Box& box, unsigned n) {
    switch (n) {
#define CASE(N) case N: return native_spectrum<N>(box)
        CASE(3); CASE(4); CASE(5); CASE(6); CASE(7); CASE(8);
#undef CASE
        default: throw std::runtime_error("PEIGEN dimension out of range");
    }
}
}
