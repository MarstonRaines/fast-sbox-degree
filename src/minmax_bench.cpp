#include "common.hpp"
#include <atomic>
int main(int argc,char** argv) {
    try {
        require(argc==6,"METHOD min|max INPUT SWAPS STEPS");
        auto method=method_of(argv[1]);auto metric=metric_of(argv[2]);
        require(metric==Metric::Minimum||metric==Metric::Maximum,"min|max required");
        auto input=read_input(argv[3]);auto swaps=read_swaps(argv[4],input.box.size());
        uint64_t steps=std::stoull(argv[5]);require(steps>0,"positive steps required");
        warm_cpu();auto begin=Clock::now();Context c(input.n,input.m);auto ready=Clock::now();
        Engine e(c,method,metric,input.box,true);auto initialized=Clock::now();uint64_t digest=0;
        for(uint64_t i=0;i<steps;++i) {
            // Force each timed call to observe the prepared input afresh.
            std::atomic_signal_fence(std::memory_order_seq_cst);
            if(metric==Metric::Minimum)e.value=e.evaluate();
            else {auto[p,q]=swaps[i%swaps.size()];e.apply(p,q,false);}
            digest=fold(digest,e.value);
        }
        auto end=Clock::now();
        std::cout<<"{\"context_ns\":"<<elapsed(begin,ready)<<",\"init_ns\":"<<elapsed(ready,initialized)
          <<",\"kernel_ns\":"<<elapsed(initialized,end)<<",\"steps\":"<<steps<<",\"final_value\":"<<e.value
          <<",\"final_box_hash\":"<<box_hash(e.box)<<",\"digest\":"<<digest<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
