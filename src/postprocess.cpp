// Reuse the frozen benchmark's input, oracle and JSON helpers unchanged.
#define main original_benchmark_main
#include "bench.cpp"
#undef main

// Standard vectorial NL: every nonzero output mask and every input mask.
// Early rejection is safe only above 48; equality is checked at the end.
bool nl104(const Box& box) {
    std::array<int,256> w;
    int maximum = 0;
    for (unsigned l = 1; l < 256; ++l) {
        for (unsigned x = 0; x < 256; ++x)
            w[x] = 1 - 2 * int(std::popcount(l & box[x]) & 1U);
        for (unsigned step = 1; step < 256; step *= 2)
            for (unsigned base = 0; base < 256; base += 2*step)
                for (unsigned j = 0; j < step; ++j) {
                    int a = w[base+j], b = w[base+step+j];
                    w[base+j] = a+b; w[base+step+j] = a-b;
                }
        for (int v : w) {
            maximum = std::max(maximum, std::abs(v));
            if (maximum > 48) return false;
        }
    }
    return maximum == 48;
}

int direct_nl(const Box& box) {
    int maximum = 0;
    for (unsigned l = 1; l < 256; ++l)
        for (unsigned a = 0; a < 256; ++a) {
            int v = 0;
            for (unsigned x = 0; x < 256; ++x)
                v += 1-2*int((std::popcount(l & box[x])+std::popcount(a & x)) & 1U);
            maximum = std::max(maximum, std::abs(v));
        }
    return 128-maximum/2;
}

using Score = std::pair<unsigned,uint64_t>;
Score score(const Engine& e) { return {hist_min(e.hist), spectrum_sum(e.hist)}; }
struct Step { unsigned p,q; uint64_t attempt; Score value; };
struct Run {
    uint64_t attempts=0, nl_calls=0, degree_rejects=0, nl_rejects=0;
    uint64_t digest=14695981039346656037ULL;
    std::vector<Step> accepted;
};

// Same first-improvement/restart skeleton as the existing search(), now with
// the predeclared lexicographic degree objective and exact NL constraint.
Run constrained_search(Engine& e, Engine* paired=nullptr) {
    Run r;
    while (score(e).first < 7) {
        bool improved=false;
        Score current=score(e);
        for (unsigned p=0; p<255 && !improved; ++p) {
            for (unsigned q=p+1; q<256; ++q) {
                Box old_box;
                std::vector<Word> old_anf, old_lanes, old_truth;
                std::vector<uint8_t> old_degrees;
                Hist old_hist=e.hist;
                if (paired) {
                    old_box=e.box; old_anf=e.component_anfs;
                    old_lanes=e.degree_lanes; old_degrees=e.degrees;
                    old_truth=paired->truth;
                }
                e.apply(p,q,true);
                Score candidate=score(e);
                if (paired) {
                    paired->apply(p,q,true);
                    require(e.hist==paired->hist && e.box==paired->box &&
                            candidate==score(*paired), "paired candidate mismatch");
                }
                ++r.attempts;
                r.digest=fold(fold(r.digest,candidate.first),candidate.second);
                bool accept=false;
                if (candidate>current) {
                    ++r.nl_calls;
                    accept=nl104(e.box);
                    if (!accept) ++r.nl_rejects;
                } else ++r.degree_rejects;
                r.digest=fold(r.digest,accept);
                if (accept) {
                    if (paired) {
                        require(direct_nl(e.box)==104, "accepted NL oracle mismatch");
                        require(from_degrees(oracle_degrees(e.c,e.box))==e.hist,
                                "accepted ANF oracle mismatch");
                    }
                    r.accepted.push_back({p,q,r.attempts,candidate});
                    improved=true;
                    break;
                }
                e.undo(p,q);
                if (paired) {
                    paired->undo(p,q);
                    require(e.box==old_box && paired->box==old_box &&
                            e.hist==old_hist && paired->hist==old_hist &&
                            score(e)==current && score(*paired)==current &&
                            e.component_anfs==old_anf && e.degree_lanes==old_lanes &&
                            e.degrees==old_degrees && paired->truth==old_truth,
                            "paired rollback state mismatch");
                }
            }
        }
        if (!improved) break;
    }
    return r;
}

int main(int argc, char** argv) {
    try {
        require(argc>=2,"check-nl AES_LUT | verify INPUT | bench METHOD INPUT");
        std::string mode=argv[1];
        if (mode=="check-nl") {
            require(argc==3,"check-nl AES_LUT");
            auto aes=read_input(argv[2]);
            require(aes.n==8 && aes.m==8,"expected 8-bit AES");
            Box identity(256); std::iota(identity.begin(),identity.end(),0);
            require(direct_nl(identity)==0 && !nl104(identity),"NL<104 rejection failed");
            require(direct_nl(aes.box)==112 && !nl104(aes.box),"NL>104 rejection failed");
            std::cout << "{\"status\":\"pass\",\"rejected_NL\":[0,112]}\n";
            return 0;
        }
        bool verify=mode=="verify";
        require((verify && argc==3) || (mode=="bench" && argc==4),"invalid command");
        auto input=read_input(argv[verify?2:3]);
        require(input.n==8 && input.m==8,"postprocessing requires an 8-bit permutation");
        Box sorted=input.box; std::sort(sorted.begin(),sorted.end());
        for (unsigned i=0; i<256; ++i) require(sorted[i]==i,"input is not a permutation");
        Method method=verify?Method::Proposed:method_of(argv[2]);
        require(method==Method::Proposed || method==Method::Peigen,"invalid evaluator");
        require(nl104(input.box),"initial NL is not exactly 104");
        auto start=Clock::now();
        Context c(8,8);
        Engine e(c,method,Metric::Spectrum,input.box,true);
        auto initial=score(e);
        std::optional<Engine> paired;
        if (verify) paired.emplace(c,Method::Peigen,Metric::Spectrum,input.box,true);
        if (verify) {
            require(e.hist==paired->hist && direct_nl(input.box)==104,
                    "initial paired/oracle mismatch");
        }
        auto run=constrained_search(e,paired?&*paired:nullptr);
        auto end=Clock::now();
        require(e.hist==from_degrees(oracle_degrees(c,e.box)),"final ANF oracle mismatch");
        require(direct_nl(e.box)==104,"final exact NL oracle mismatch");
        require(run.degree_rejects && run.nl_rejects && !run.accepted.empty(),
                "missing acceptance/rejection branch coverage");
        std::cout << "{\"verified\":true,\"paired\":" << (verify?"true":"false")
                  << ",\"method\":" << std::quoted(verify?"paired":argv[2])
                  << ",\"total_ns\":" << ns(start,end)
                  << ",\"attempts\":" << run.attempts
                  << ",\"nl_calls\":" << run.nl_calls
                  << ",\"degree_rejects\":" << run.degree_rejects
                  << ",\"nl_rejects\":" << run.nl_rejects
                  << ",\"candidate_digest\":" << run.digest
                  << ",\"initial_min\":" << initial.first
                  << ",\"initial_sum\":" << initial.second
                  << ",\"final_min\":" << score(e).first
                  << ",\"final_sum\":" << score(e).second
                  << ",\"accepted\":[";
        for (size_t i=0; i<run.accepted.size(); ++i) {
            if (i) std::cout << ',';
            auto s=run.accepted[i];
            std::cout << '[' << s.p << ',' << s.q << ',' << s.attempt
                      << ',' << s.value.first << ',' << s.value.second << ']';
        }
        std::cout << "],\"final_histogram\":"; json_array(std::cout,e.hist);
        std::cout << ",\"final_box\":"; json_array(std::cout,e.box);
        std::cout << "}\n";
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n'; return 1;
    }
}
