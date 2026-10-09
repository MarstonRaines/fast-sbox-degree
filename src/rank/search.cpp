#include "search_constraints.hpp"

int direct_nl(const Box& box) {
    int maximum=0;
    for(unsigned l=1;l<256;++l)for(unsigned a=0;a<256;++a){int v=0;
        for(unsigned x=0;x<256;++x)v+=1-2*int((std::popcount(l&box[x])+std::popcount(a&x))&1U);
        maximum=std::max(maximum,std::abs(v));}
    return 128-maximum/2;
}
unsigned direct_du(const Box& box) {
    unsigned maximum=0;
    for(unsigned a=1;a<256;++a){std::array<unsigned,256> counts{};
        for(unsigned x=0;x<256;++x)maximum=std::max(maximum,++counts[box[x]^box[x^a]]);}
    return maximum;
}
using Score=std::pair<unsigned,uint64_t>;
template<class E> Score score(const E& e) {
    // The PEIGEN group scores minimum degree only; no spectrum tie-break.
    if constexpr(std::is_same_v<E,Engine>)return {unsigned(e.value),0};
    else return {hist_min(e.hist),spectrum_sum(e.hist)};
}
struct StaticSpectrumEngine : RankEngine {
    StaticSpectrumEngine(const Context& c,const Box& input):RankEngine(c,input,true){}
    void undo(unsigned p,unsigned q) {std::swap(box[p],box[q]);hist=old_hist;value=old_value;}
};
struct Step {unsigned p,q;uint64_t attempt;Score value;int nl;unsigned du;};
struct Run {
    uint64_t attempts=0,nl_calls=0,du_calls=0,degree_rejects=0,nl_rejects=0,du_rejects=0;
    uint64_t digest=14695981039346656037ULL;
    uint64_t degree_ns=0,nl_ns=0,du_ns=0,rollback_ns=0;
    std::vector<Step> accepted;
};

template<bool Diagnostic,bool Verify,class E> Run search(E& e,SearchConstraints& constraints) {
    Run r;std::unique_ptr<Engine> paired;
    const int target_nl=constraints.nonlinearity();const unsigned target_du=constraints.uniformity();
    if constexpr(Verify)paired=std::make_unique<Engine>(e.c,Method::Peigen,Metric::Spectrum,e.box,true);
    while(score(e).first<7) {
        bool improved=false;Score current=score(e);
        for(unsigned p=0;p<255&&!improved;++p)for(unsigned q=p+1;q<256;++q) {
            const unsigned old_p=e.box[p],old_q=e.box[q];Clock::time_point begin;
            if constexpr(Diagnostic)begin=Clock::now();
            e.apply(p,q,true);
            if constexpr(Diagnostic)r.degree_ns+=elapsed(begin,Clock::now());
            Score candidate=score(e);
            if constexpr(Verify) {
                paired->apply(p,q,true);
                require(candidate.first==hist_min(paired->hist)&&e.box==paired->box,"candidate degree/state mismatch");
                if constexpr(!std::is_same_v<E,Engine>)require(e.hist==paired->hist,"candidate spectrum mismatch");
            }
            ++r.attempts;r.digest=fold(fold(r.digest,candidate.first),candidate.second);
            bool nl_changed=false,du_changed=false,accept=false;
            if(candidate>current) {
                ++r.nl_calls;nl_changed=true;
                if constexpr(Diagnostic)begin=Clock::now();
                constraints.menyachikhin_walsh_update(p,q,old_p,old_q);
                const bool nl_ok=constraints.nonlinearity()==target_nl;
                if constexpr(Diagnostic)r.nl_ns+=elapsed(begin,Clock::now());
                if constexpr(Verify)constraints.verify(e.box,true,false);
                if(nl_ok) {
                    ++r.du_calls;du_changed=true;
                    if constexpr(Diagnostic)begin=Clock::now();
                    constraints.menyachikhin_ddt_update(e.box,p,q);
                    accept=constraints.uniformity()==target_du;
                    if constexpr(Diagnostic)r.du_ns+=elapsed(begin,Clock::now());
                    if constexpr(Verify)constraints.verify(e.box,false,true);
                    if(!accept)++r.du_rejects;
                } else ++r.nl_rejects;
            } else ++r.degree_rejects;
            r.digest=fold(fold(fold(r.digest,nl_changed),du_changed),accept);
            if(accept) {
                r.accepted.push_back({p,q,r.attempts,candidate,constraints.nonlinearity(),constraints.uniformity()});
                improved=true;break;
            }
            if constexpr(Diagnostic)begin=Clock::now();
            if(du_changed)constraints.menyachikhin_ddt_update(e.box,p,q,true);
            if(nl_changed)constraints.menyachikhin_walsh_update(p,q,old_q,old_p);
            e.undo(p,q);
            if constexpr(Diagnostic)r.rollback_ns+=elapsed(begin,Clock::now());
            if constexpr(Verify) {
                paired->undo(p,q);require(score(e)==current&&e.box==paired->box,"degree rollback mismatch");
                if constexpr(!std::is_same_v<E,Engine>)require(e.hist==paired->hist,"spectrum rollback mismatch");
                if(nl_changed)constraints.verify(e.box);
            }
        }
        if(!improved)break;
    }
    return r;
}

template<bool Diagnostic,bool Verify,class E> void run(const Input& input,const std::string& method) {
    require(input.n==8&&input.m==8,"8-bit input required");
    Box sorted=input.box;std::sort(sorted.begin(),sorted.end());for(unsigned i=0;i<256;++i)require(sorted[i]==i,"input not a permutation");
    warm_cpu();auto begin=Clock::now();Context c(8,8);auto context_end=Clock::now();
    std::unique_ptr<E> e;
    if constexpr(std::is_base_of_v<RankEngine,E>)e=std::make_unique<E>(c,input.box);
    else e=std::make_unique<E>(c,Method::Peigen,Metric::Minimum,input.box,true);
    auto degree_init_end=Clock::now();SearchConstraints constraints(8,input.box);auto init_end=Clock::now();
    require(constraints.nonlinearity()==104&&constraints.uniformity()==8,"initial NL104/DU8 required");
    auto initial=score(*e);auto r=search<Diagnostic,Verify>(*e,constraints);auto end=Clock::now();
    // Independent accepted-state recomputation is outside timing. Its spectrum
    // is used for reporting only, never for the minimum-only search decisions.
    Box verified=input.box;auto initial_hist=independent_histogram(c,verified);
    require(initial.first==hist_min(initial_hist)&&direct_nl(verified)==104&&direct_du(verified)==8,"initial independent properties");
    std::vector<Hist> accepted_hist;
    for(const auto& a:r.accepted) {
        std::swap(verified[a.p],verified[a.q]);auto hist=independent_histogram(c,verified);
        require(a.value.first==hist_min(hist)&&direct_nl(verified)==a.nl&&a.nl==104&&direct_du(verified)==a.du&&a.du==8,"accepted independent properties");
        if constexpr(!std::is_same_v<E,Engine>)require(a.value.second==spectrum_sum(hist),"accepted spectrum sum");
        accepted_hist.push_back(hist);
    }
    const auto final_hist=independent_histogram(c,e->box);
    require(verified==e->box&&hist_min(final_hist)==score(*e).first,"final independent state");
    if constexpr(!std::is_same_v<E,Engine>)require(e->hist==final_hist,"final spectrum");
    constraints.verify(e->box);const bool reached=score(*e).first==7;
    std::cout<<"{\"protocol\":\"fixed-nl-du\",\"method\":"<<std::quoted(method)
      <<",\"objective\":"<<std::quoted(std::is_same_v<E,Engine>?"minimum-only":"minimum-then-spectrum-sum")
      <<",\"diagnostic\":"<<(Diagnostic?"true":"false")<<",\"all_candidates_paired\":"<<(Verify?"true":"false")
      <<",\"context_ns\":"<<elapsed(begin,context_end)<<",\"degree_init_ns\":"<<elapsed(context_end,degree_init_end)
      <<",\"constraints_init_ns\":"<<elapsed(degree_init_end,init_end)<<",\"init_ns\":"<<elapsed(context_end,init_end)
      <<",\"kernel_ns\":"<<elapsed(init_end,end)<<",\"total_ns\":"<<elapsed(begin,end)
      <<",\"attempts\":"<<r.attempts<<",\"nl_calls\":"<<r.nl_calls<<",\"du_calls\":"<<r.du_calls
      <<",\"degree_rejects\":"<<r.degree_rejects<<",\"nl_rejects\":"<<r.nl_rejects<<",\"du_rejects\":"<<r.du_rejects
      <<",\"candidate_digest\":"<<r.digest<<",\"initial_min\":"<<hist_min(initial_hist)<<",\"initial_sum\":"<<spectrum_sum(initial_hist)
      <<",\"final_min\":"<<hist_min(final_hist)<<",\"final_sum\":"<<spectrum_sum(final_hist)
      <<",\"initial_nl\":104,\"initial_du\":8,\"final_nl\":"<<constraints.nonlinearity()<<",\"final_du\":"<<constraints.uniformity()
      <<",\"target_reached\":"<<(reached?"true":"false")<<",\"termination\":"<<std::quoted(reached?"target-reached":"no-strict-improvement")
      <<",\"independently_verified_states\":"<<r.accepted.size()+1<<",\"accepted\":[";
    for(size_t i=0;i<r.accepted.size();++i) {
        if(i)std::cout<<',';
        const auto& a=r.accepted[i];
        std::cout<<"{\"p\":"<<a.p<<",\"q\":"<<a.q<<",\"attempt\":"<<a.attempt<<",\"minimum\":"<<a.value.first
          <<",\"spectrum_sum\":"<<spectrum_sum(accepted_hist[i])<<",\"nonlinearity\":"<<a.nl<<",\"differential_uniformity\":"<<a.du<<'}';
    }
    std::cout<<"],\"final_histogram\":";json_array(final_hist);std::cout<<",\"final_box\":";json_array(e->box);
    if constexpr(Diagnostic)std::cout<<",\"degree_ns\":"<<r.degree_ns<<",\"nl_ns\":"<<r.nl_ns<<",\"du_ns\":"<<r.du_ns
       <<",\"rollback_ns\":"<<r.rollback_ns<<",\"other_ns\":"<<elapsed(init_end,end)-r.degree_ns-r.nl_ns-r.du_ns-r.rollback_ns;
    std::cout<<"}\n";
}
template<bool Diagnostic,bool Verify> void dispatch(const Input& input,const std::string& method) {
    if(method=="rank-shared")run<Diagnostic,Verify,RankEngine>(input,method);
    else if(method=="rank-static")run<Diagnostic,Verify,StaticSpectrumEngine>(input,method);
    else {require(method=="peigen","method");run<Diagnostic,Verify,Engine>(input,method);}
}
int main(int argc,char** argv) {
    try {
        require(argc==4,"bench|diagnostic|verify METHOD INPUT");std::string mode=argv[1],method=argv[2];auto input=read_input(argv[3]);
        if(mode=="verify")dispatch<false,true>(input,method);
        else if(mode=="diagnostic")dispatch<true,false>(input,method);
        else {require(mode=="bench","mode");dispatch<false,false>(input,method);}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
