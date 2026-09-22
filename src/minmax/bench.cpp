#include "degree.hpp"
#include "reference.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <type_traits>

using namespace coam;
using Clock = std::chrono::steady_clock;
using Swap = std::pair<unsigned, unsigned>;
uint64_t ns(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b-a).count();
}
void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
uint64_t fold(uint64_t h, uint64_t x) { return (h ^ x) * 1099511628211ULL; }
uint64_t box_hash(const Box& s) {
    uint64_t h = 14695981039346656037ULL;
    for (auto x : s) h = fold(h, x);
    return h;
}
template<class T> void json_array(std::ostream& out, const T& a) {
    out << '[';
    bool first = true;
    for (auto x : a) { if (!first) out << ','; first = false; out << +x; }
    out << ']';
}
struct Input { unsigned n, m; Box box; };
Input read_input(const std::string& path) {
    std::ifstream in(path);
    Input x{};
    require(bool(in >> x.n >> x.m), "cannot read input: " + path);
    require(x.n >= 1 && x.n <= 20 && x.m >= 1 && x.m <= 20, "invalid dimensions");
    x.box.resize(1U << x.n);
    for (auto& v : x.box) require(bool(in >> v) && v < (1U << x.m), "invalid LUT entry");
    std::string extra;
    require(!(in >> extra), "unexpected extra LUT data");
    return x;
}
std::vector<Swap> read_swaps(const std::string& path, unsigned length) {
    std::ifstream in(path);
    require(bool(in), "cannot read swaps: " + path);
    std::vector<Swap> pairs;
    unsigned p, q;
    while (in >> p) {
        require(bool(in >> q) && p < length && q < length && p != q, "invalid swap");
        pairs.emplace_back(p, q);
    }
    require(in.eof() && !pairs.empty(), "invalid or empty swap file");
    return pairs;
}

// Independent byte-array truth-table transform, used only for correctness.
std::vector<unsigned> oracle_degrees(const Context& c, const Box& box) {
    std::vector<unsigned> result(c.components);
    std::vector<uint8_t> anf(c.length);
    for (unsigned l = 1; l < c.components; ++l) {
        for (unsigned x = 0; x < c.length; ++x) anf[x] = std::popcount(l & box[x]) & 1U;
        for (unsigned step = 1; step < c.length; step *= 2)
            for (unsigned base = 0; base < c.length; base += 2 * step)
                for (unsigned j = 0; j < step; ++j) anf[base + step + j] ^= anf[base + j];
        for (unsigned u = 0; u < c.length; ++u)
            if (anf[u]) result[l] = std::max(result[l], unsigned(std::popcount(u)));
    }
    return result;
}
Hist from_degrees(const std::vector<unsigned>& deg) {
    Hist h{};
    for (size_t l = 1; l < deg.size(); ++l) ++h[deg[l]];
    return h;
}

void verify_box(const Context& c, const Box& box, bool use_oracle) {
    Engine bit(c, Method::Bitwise, Metric::Spectrum, box);
    Engine prop(c, Method::Proposed, Metric::Spectrum, box);
    Engine mini(c, Method::Proposed, Metric::Minimum, box);
    Engine maxi(c, Method::Proposed, Metric::Maximum, box);
    require(bit.hist == prop.hist, "spectrum disagreement");
    require(mini.value == hist_min(bit.hist) && maxi.value == hist_max(bit.hist), "extremum disagreement");
    std::vector<unsigned> reference;
    if (use_oracle) {
        reference = oracle_degrees(c, box);
        require(bit.hist == from_degrees(reference), "independent oracle disagreement");
        std::vector<unsigned> expected;
        for (unsigned l = 1; l < c.components; ++l)
            if (reference[l] == mini.value) expected.push_back(l);
        require(mini.minimum_lambdas() == expected, "minimum lambda enumeration disagreement");
    }
    if (c.n == c.m && c.n >= 3 && c.n <= 8) {
        require(peigen_native_spectrum(box, c.n) == bit.hist, "native PEIGEN spectrum disagreement");
        for (Metric metric : {Metric::Minimum, Metric::Maximum, Metric::Spectrum}) {
            Hist h{};
            auto d = peigen_evaluate(box, c.n, metric, h);
            require(d == (metric == Metric::Minimum ? mini.value : metric == Metric::Maximum ? maxi.value : bit.value),
                    "PEIGEN numerical core disagreement");
            if (metric == Metric::Spectrum) require(h == bit.hist, "PEIGEN histogram disagreement");
        }
    }
}

void verify_updates(const Context& c, const Box& initial, const std::vector<Swap>& swaps) {
    Engine spectrum(c, Method::Proposed, Metric::Spectrum, initial);
    Engine maximum(c, Method::Proposed, Metric::Maximum, initial);
    Engine packed_bit(c, Method::Bitwise, Metric::Spectrum, initial, true);
    std::optional<Engine> packed_peigen;
    if (c.n == c.m && c.n >= 3 && c.n <= 8)
        packed_peigen.emplace(c,Method::Peigen,Metric::Spectrum,initial,true);
    Box current = initial;
    uint64_t index = 0;
    for (auto [p,q] : swaps) {
        Hist previous_hist = spectrum.hist;
        auto previous_degrees = spectrum.degrees;
        auto previous_max_degrees = maximum.degrees;
        // Byte-for-byte state checks are bounded to small inputs; all sizes check
        // full histograms after both accepted updates and rollbacks.
        std::vector<Word> previous_anfs;
        if (c.n <= 8) previous_anfs = spectrum.component_anfs;
        auto previous_coordinates = maximum.coord;
        std::swap(current[p], current[q]);
        spectrum.apply(p, q, true);
        maximum.apply(p, q, true);
        packed_bit.apply(p,q,true);
        if (packed_peigen) packed_peigen->apply(p,q,true);
        Engine bit(c, Method::Bitwise, Metric::Spectrum, current);
        Engine mini(c, Method::Proposed, Metric::Minimum, current);
        require(spectrum.hist == bit.hist, "update histogram mismatch");
        require(packed_bit.hist == bit.hist, "prepared input update mismatch");
        if (packed_peigen) require(packed_peigen->hist == bit.hist, "prepared PEIGEN update mismatch");
        require(maximum.value == hist_max(bit.hist) && mini.value == hist_min(bit.hist), "update extremum mismatch");
        // The same sequence exercises commitment and rejection.
        if ((index++ % 3) != 0) {
            std::swap(current[p], current[q]);
            spectrum.undo(p,q);
            maximum.undo(p,q);
            packed_bit.undo(p,q);
            if (packed_peigen) packed_peigen->undo(p,q);
            require(packed_bit.hist == previous_hist,
                    "prepared input rollback mismatch");
            if (packed_peigen) require(packed_peigen->hist == previous_hist, "prepared PEIGEN rollback mismatch");
            require(spectrum.box == current && maximum.box == current, "rollback LUT mismatch");
            require(spectrum.hist == previous_hist && spectrum.degrees == previous_degrees,
                    "rollback degree state mismatch");
            require(maximum.coord == previous_coordinates && maximum.degrees == previous_max_degrees,
                    "rollback coordinate state mismatch");
            if (c.n <= 8) require(spectrum.component_anfs == previous_anfs, "rollback ANF mismatch");
        }
    }
}

void verify_reference(const Context& c, const Box& initial, const std::vector<Swap>& swaps) {
    for (Metric metric : {Metric::Minimum,Metric::Maximum,Metric::Spectrum}) {
        ReferenceEngine a(c,Method::Reference,metric,initial);
        ReferenceEngine p(c,Method::ReferenceProposed,metric,initial);
        require(a.value == p.value,"reference initial scalar mismatch");
        if (metric == Metric::Spectrum) require(a.hist == p.hist,"reference initial histogram mismatch");
        for (size_t i = 0; i < swaps.size(); ++i) {
            auto [x,y] = swaps[i];
            Box previous = p.box;
            Hist previous_hist = p.hist;
            auto previous_degrees = p.degrees;
            auto previous_coord = p.coord;
            std::vector<int32_t> previous_anfs;
            if (c.n <= 8) previous_anfs = p.component_anfs;
            require(a.apply(x,y,true) == p.apply(x,y,true),"reference update scalar mismatch");
            Engine oracle(c,Method::Bitwise,Metric::Spectrum,p.box);
            require(p.value == (metric == Metric::Minimum ? hist_min(oracle.hist) :
                    metric == Metric::Maximum ? hist_max(oracle.hist) : oracle.value),
                    "reference recomputation mismatch");
            if (metric == Metric::Spectrum)
                require(a.hist == p.hist && p.hist == oracle.hist,"reference update full histogram mismatch");
            if (metric != Metric::Minimum && p.delta_output) {
                uint64_t expected = (1ULL << (c.n-std::popcount(x))) +
                    (1ULL << (c.n-std::popcount(y))) - (2ULL << (c.n-std::popcount(x|y)));
                require(p.delta.size() == expected,"symmetric difference size mismatch");
                require(p.top_delta.size() == unsigned(std::popcount(x^y)),"top difference size mismatch");
            }
            if (i % 3 != 0) {
                a.undo(x,y); p.undo(x,y);
                require(a.box == previous && p.box == previous,"reference rollback LUT mismatch");
                require(p.degrees == previous_degrees && p.hist == previous_hist,"reference rollback degrees mismatch");
                if (metric != Metric::Minimum) require(p.coord == previous_coord,"reference rollback coordinates mismatch");
                if (c.n <= 8) require(p.component_anfs == previous_anfs,"reference rollback ANFs mismatch");
                require(a.evaluate() == p.value,"reference rollback value mismatch");
            }
        }
    }
}

void self_check() {
    uint64_t cases = 0, swaps = 0;
    Context small(2,2);
    std::vector<Swap> all_pairs;
    for (unsigned p = 0; p < 4; ++p) for (unsigned q = p+1; q < 4; ++q) all_pairs.emplace_back(p,q);
    for (unsigned encoded = 0; encoded < 256; ++encoded) {
        Box box(4);
        for (unsigned i = 0; i < 4; ++i) box[i] = (encoded >> (2*i)) & 3U;
        verify_box(small, box, true);
        verify_updates(small, box, all_pairs);
        verify_reference(small,box,all_pairs);
        ++cases; swaps += all_pairs.size();
    }
    std::mt19937 gen(20260917);
    for (unsigned n = 1; n <= 9; ++n) {
        for (unsigned m : {n, std::min(10U, n+1)}) {
            Context c(n,m);
            for (unsigned type = 0; type < 5; ++type) {
                Box box(c.length);
                for (unsigned x = 0; x < c.length; ++x) {
                    if (type == 0) box[x] = 0;
                    if (type == 1) box[x] = c.components-1;
                    if (type == 2) box[x] = ((x & 1U) ? 3U : 0U) & (c.components-1);
                    if (type == 3) box[x] = (4U | ((x & 1U) ? 3U : 0U)) & (c.components-1);
                    if (type == 4) box[x] = gen() & (c.components-1);
                }
                verify_box(c, box, true);
                std::vector<Swap> pairs;
                for (unsigned j = 0; j < 12; ++j) {
                    unsigned p = gen() % c.length, q = gen() % (c.length-1);
                    if (q >= p) ++q;
                    pairs.emplace_back(p,q);
                }
                verify_updates(c,box,pairs);
                verify_reference(c,box,pairs);
                ++cases; swaps += pairs.size();
            }
        }
    }
#ifdef COAM_COUNT
    require(counts.top_guard_kept && counts.top_guard_fallback && counts.degree_drops,
            "packed guard branch coverage missing");
    require(counts.scalar_guard_kept && counts.scalar_guard_fallback && counts.scalar_degree_drops,
            "reference guard branch coverage missing");
#endif
    std::cout << "{\"status\":\"pass\",\"static_cases\":" << cases
              << ",\"swap_cases_per_implementation_group\":" << swaps << ",\"exhaustive_2x2_maps\":256,\"counts\":";
    print_counts(std::cout);
    std::cout << "}\n";
}

struct Accepted { unsigned p,q; uint64_t value; };
struct SearchResult {
    uint64_t attempts = 0, digest = 14695981039346656037ULL;
    std::vector<Accepted> accepted;
};
template<class E> SearchResult search(E& e, E* paired = nullptr) {
    SearchResult r;
    while (true) {
        bool improved = false;
        uint64_t current = e.value;
        for (unsigned p = 0; p < e.c.length && !improved; ++p) {
            for (unsigned q = p+1; q < e.c.length; ++q) {
                uint64_t candidate = e.apply(p,q,true);
                if (paired) {
                    require(paired->apply(p,q,true) == candidate, "paired search candidate score mismatch");
                    if (e.metric == Metric::Spectrum)
                        require(e.hist == paired->hist, "paired search candidate histogram mismatch");
                }
                ++r.attempts;
                r.digest = fold(r.digest, candidate);
                if (candidate > current) {
                    r.accepted.push_back({p,q,candidate});
                    improved = true;
                    break;
                }
                e.undo(p,q);
                if (paired) paired->undo(p,q);
            }
        }
        if (!improved) break;
    }
    return r;
}

template<class E> void benchmark(int argc, char** argv) {
    require(argc == 9, "bench METHOD METRIC MODE INPUT SWAPS LOOPS TRACE_JSON");
    std::string method_name = argv[2], metric_name = argv[3], mode = argv[4];
    Method method = method_of(method_name);
    Metric metric = metric_of(metric_name);
    auto input = read_input(argv[5]);
    require((input.n <= 16 && input.m <= 16) ||
            (method == Method::Proposed || method == Method::Bitwise),
            "above 16 bits only proposed/bitwise benchmarks are supported");
    auto loops = std::stoull(argv[7]);
    require(loops > 0, "loops must be positive");
    require(mode == "static" || mode == "updates" || mode == "search", "unknown mode");
    require((mode != "static" || metric == Metric::Minimum), "static benchmark is minimum evaluation");
    require(mode != "search" || loops == 1, "search runs one complete trajectory");
    std::vector<Swap> pairs;
    if (mode == "updates") pairs = read_swaps(argv[6], input.box.size());
    auto a = Clock::now();
    Context c(input.n,input.m);
    auto b = Clock::now();
    E e(c, method, metric, input.box, mode != "static");
    counts = {};
    auto begin = Clock::now();
    uint64_t initial = e.value, steps = 0, digest = 14695981039346656037ULL;
    SearchResult sr;
    if (mode == "search") {
        sr = search(e);
        steps = sr.attempts;
        digest = sr.digest;
    } else if (mode == "static") {
        for (uint64_t i = 0; i < loops; ++i) {
            asm volatile("" : : "g"(e.box.data()) : "memory");
            e.value = e.evaluate();
            digest = fold(digest, e.value);
            ++steps;
        }
    } else {
        for (uint64_t i = 0; i < loops; ++i) for (auto [p,q] : pairs) {
            e.apply(p,q,false);
            digest = fold(digest,e.value);
            ++steps;
        }
    }
    auto end = Clock::now();
    Counts kernel_counts = counts;
    // All checks, formatting, and file I/O are outside the timed interval.
    Engine final_check(c, Method::Bitwise, Metric::Spectrum, e.box);
    require(e.value == (metric == Metric::Minimum ? hist_min(final_check.hist) :
                       metric == Metric::Maximum ? hist_max(final_check.hist) : final_check.value),
            "final recomputation check failed");
    if (metric == Metric::Spectrum) require(e.hist == final_check.hist, "final full histogram check failed");
    counts = kernel_counts;
    std::ostringstream record;
    record << "{\"method\":" << std::quoted(method_name) << ",\"metric\":" << std::quoted(metric_name)
           << ",\"mode\":" << std::quoted(mode) << ",\"n\":" << c.n << ",\"m\":" << c.m
           << ",\"loops\":" << loops << ",\"steps\":" << steps << ",\"accepted\":" << sr.accepted.size()
           << ",\"context_ns\":" << ns(a,b) << ",\"init_ns\":" << ns(b,begin)
           << ",\"kernel_ns\":" << ns(begin,end) << ",\"total_ns\":" << ns(a,end)
           << ",\"initial_value\":" << initial << ",\"final_value\":" << e.value
           << ",\"candidate_digest\":" << digest << ",\"final_box_hash\":" << box_hash(e.box)
           << ",\"final_histogram\":";
    json_array(record, final_check.hist);
    record << ",\"verified\":true";
    if constexpr (std::is_same_v<E,ReferenceEngine>)
        record << ",\"input_conversion_ns\":" << e.conversion_ns
               << ",\"state_init_ns\":" << ns(b,begin)-e.conversion_ns;
#ifdef COAM_COUNT
    record << ",\"instrumented\":true,\"counts\":";
    print_counts(record);
#else
    record << ",\"instrumented\":false";
#endif
    std::ofstream trace(argv[8]);
    require(bool(trace), "cannot create trace file");
    trace << record.str() << ",\"accepted_swaps\":[";
    for (size_t i = 0; i < sr.accepted.size(); ++i) {
        if (i) trace << ',';
        const auto& s = sr.accepted[i];
        trace << '[' << s.p << ',' << s.q << ',' << s.value << ']';
    }
    trace << "],\"final_box\":";
    json_array(trace,e.box);
    trace << "}\n";
    require(bool(trace), "failed writing trace");
    std::cout << record.str() << "}\n";
}

#ifndef COAM_NO_BENCH_MAIN
int main(int argc, char** argv) {
    try {
        require(argc >= 2, "commands: check, verify, verify-search, bench");
        std::string command = argv[1];
        if (command == "check") self_check();
        else if (command == "verify-reference") {
            require(argc == 4,"verify-reference INPUT SWAPS");
            auto x = read_input(argv[2]);
            Context c(x.n,x.m);
            auto swaps = read_swaps(argv[3],c.length);
            verify_reference(c,x.box,swaps);
            std::cout << "{\"status\":\"pass\",\"n\":" << c.n << ",\"swaps\":" << swaps.size() << "}\n";
        } else if (command == "verify") {
            require(argc == 4, "verify INPUT SWAPS");
            auto x = read_input(argv[2]);
            Context c(x.n,x.m);
            verify_box(c,x.box,c.n <= 8);
            auto swaps = read_swaps(argv[3], c.length);
            verify_updates(c,x.box,swaps);
            std::cout << "{\"status\":\"pass\",\"n\":" << c.n << ",\"swaps\":" << swaps.size() << "}\n";
        } else if (command == "verify-search") {
            require(argc == 4, "verify-search INPUT METRIC");
            auto x = read_input(argv[2]);
            require(x.n == 8 && x.m == 8, "search protocol requires 8x8");
            Context c(x.n,x.m);
            Metric metric = metric_of(argv[3]);
            Engine p(c,Method::Proposed,metric,x.box,true), b(c,Method::Peigen,metric,x.box,true);
            require(p.value == b.value, "initial search mismatch");
            auto s = search(p,&b);
            require(p.box == b.box && p.value == b.value, "final search mismatch");
            verify_box(c,p.box,true);
            std::cout << "{\"status\":\"pass\",\"attempts\":" << s.attempts
                      << ",\"accepted\":" << s.accepted.size() << ",\"digest\":" << s.digest << "}\n";
        } else if (command == "bench") {
            require(argc >= 3,"bench requires method");
            Method method = method_of(argv[2]);
            if (method == Method::Reference || method == Method::ReferenceProposed)
                benchmark<ReferenceEngine>(argc,argv);
            else benchmark<Engine>(argc,argv);
        }
        else throw std::runtime_error("unknown command");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
#endif
