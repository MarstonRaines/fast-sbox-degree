// Verification only: reuse input parsing, without the benchmark entry point.
#define COAM_NO_BENCH_MAIN
#include "minmax/bench.cpp"

uint64_t checked_states=0, checked_coefficients=0;

void check_state(Engine& mini, Engine& maximum) {
    const Context& c=mini.c;
    // Independent unpacked vector-valued Mobius transform.
    Box coefficients=mini.box;
    for (unsigned bit=1; bit<c.length; bit*=2)
        for (unsigned u=0; u<c.length; ++u)
            if (u & bit) coefficients[u]^=coefficients[u^bit];
    std::vector<std::pair<unsigned,unsigned>> terms;
    for (unsigned u:c.order) if (coefficients[u]) terms.emplace_back(coefficients[u],c.weights[u]);
    Hist hist{};
    std::vector<unsigned> degrees(c.components), lambdas;
    for (unsigned l=1; l<c.components; ++l) {
        for (auto [v,d]:terms) if (std::popcount(l & v) & 1U) { degrees[l]=d; break; }
        ++hist[degrees[l]];
    }
    require(mini.value==hist_min(hist) && maximum.value==hist_max(hist), "independent extremum mismatch");
    for (unsigned l=1; l<c.components; ++l) if (degrees[l]==mini.value) lambdas.push_back(l);
    require(mini.minimum_lambdas()==lambdas, "lambda enumeration mismatch");
    for (unsigned u=0; u<c.length; ++u) {
        unsigned packed=c.m<=16 ? unsigned(mini.output_anfs[u/4]>>(16*(u%4))) & 65535U
                                : unsigned(mini.output_anfs[u/2]>>(32*(u%2)));
        require(packed==coefficients[u], "packed output ANF mismatch");
        for (unsigned b=0; b<c.m; ++b)
            require(((maximum.coord[b*c.words+u/64]>>(u%64)) & 1U)==((coefficients[u]>>b)&1U),
                    "maximum coordinate ANF mismatch");
        ++checked_coefficients;
    }
    Engine static_min(c,Method::Proposed,Metric::Minimum,mini.box,false);
    require(static_min.value==mini.value, "prepared/static minimum mismatch");
    ++checked_states;
}

void verify_case(unsigned n,unsigned m,const Box& box,const std::vector<Swap>& pairs) {
    Context c(n,m);
    Engine mini(c,Method::Proposed,Metric::Minimum,box,true);
    Engine maximum(c,Method::Proposed,Metric::Maximum,box,true);
    check_state(mini,maximum);
    for (size_t i=0;i<pairs.size();++i) {
        auto [p,q]=pairs[i];
        auto old_box=mini.box;
        auto old_output=mini.output_truth;
        auto old_coord=maximum.coord;
        auto old_degrees=maximum.degrees;
        auto old_value=mini.value;
        mini.apply(p,q,true); maximum.apply(p,q,true);
        check_state(mini,maximum);
        if (i%3!=0) {
            mini.undo(p,q); maximum.undo(p,q);
            require(mini.box==old_box && maximum.box==old_box && mini.output_truth==old_output &&
                    maximum.coord==old_coord && maximum.degrees==old_degrees && mini.value==old_value,
                    "rollback state mismatch");
            // Undo restores the LUT and prepared truth; evaluation refreshes temporary ANFs.
            mini.value=mini.evaluate();
            check_state(mini,maximum);
        }
    }
}

int main(int argc,char** argv) {
    try {
        require(argc>=2,"verify_minmax --edges or INPUT SWAPS");
        if (std::string(argv[1])=="--edges") {
            for (auto [n,m]:std::array<Swap,7>{{{1,17},{3,19},{6,16},{6,17},{6,20},{17,3},{19,3}}})
                for (unsigned type=0;type<4;++type) {
                    Box box(1U<<n);
                    for (unsigned x=0;x<box.size();++x)
                        box[x]=type==0 ? 0 : type==1 ? (1U<<m)-1 :
                            type==2 ? (x & 1U)*((1U<<m)-1) : (x*17U+3U)&((1U<<m)-1);
                    std::vector<Swap> pairs{{0,1},{0,unsigned(box.size()-1)}};
                    verify_case(n,m,box,pairs);
                }
        } else {
            require(argc==3,"verify_minmax INPUT SWAPS");
            auto x=read_input(argv[1]);
            auto pairs=read_swaps(argv[2],x.box.size());
            verify_case(x.n,x.m,x.box,pairs);
        }
        std::cout << "{\"status\":\"pass\",\"checked_states\":" << checked_states
                  << ",\"checked_coefficients\":" << checked_coefficients << "}\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
