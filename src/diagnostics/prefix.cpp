#include "common.hpp"

// Separate prefix-cost experiment: checkpoint overhead never enters the main suite.
template<class E> void prefix_run(const Input& input,const std::string& method,
        const std::vector<Swap>& swaps,unsigned maximum) {
    std::vector<unsigned> points;
    for(unsigned k:{1U,16U,128U,512U})if(k<=maximum)points.push_back(k);
    require(!points.empty()&&points.back()==maximum,"maximum must be 1,16,128,512");
    std::vector<uint64_t> times(points.size()),digests(points.size());
    std::vector<Hist> histograms(points.size());
    warm_cpu();auto begin=Clock::now();Context c(input.n,input.m);auto context_end=Clock::now();
    std::unique_ptr<E> e;
    if constexpr(std::is_same_v<E,RankEngine>)e=std::make_unique<E>(c,input.box,method=="rank-static");
    else e=std::make_unique<E>(c,method_of(method),Metric::Spectrum,input.box,true);
    auto init_end=Clock::now();uint64_t digest=14695981039346656037ULL;size_t point=0;
    for(unsigned i=1;i<=maximum;++i) {
        auto[p,q]=swaps[(i-1)%swaps.size()];digest=fold(digest,e->apply(p,q,false));
        if(i==points[point]) {
            times[point]=elapsed(init_end,Clock::now());digests[point]=digest;histograms[point]=e->hist;++point;
        }
    }
    std::cout<<"{\"method\":"<<std::quoted(method)<<",\"n\":"<<input.n<<",\"m\":"<<input.m
       <<",\"context_ns\":"<<elapsed(begin,context_end)<<",\"init_ns\":"<<elapsed(context_end,init_end)
       <<",\"peak_rss_bytes\":"<<peak_rss_bytes()<<",\"context_allocated_bytes\":"<<context_storage(c)
       <<",\"state_allocated_bytes\":"<<state_storage(*e)<<",\"checkpoints\":[";
    for(size_t i=0;i<points.size();++i){if(i)std::cout<<',';
        std::cout<<"{\"steps\":"<<points[i]<<",\"kernel_ns\":"<<times[i]
           <<",\"total_ns\":"<<elapsed(begin,init_end)+times[i]<<",\"digest\":"<<digests[i]<<",\"histogram\":";
        json_array(histograms[i]);std::cout<<'}';}
    std::cout<<"]}\n";
}
int main(int argc,char** argv) {
    try{
        require(argc==5,"METHOD INPUT SWAPS MAXIMUM");std::string method=argv[1];auto input=read_input(argv[2]);
        auto swaps=read_swaps(argv[3],input.box.size());unsigned maximum=std::stoul(argv[4]);
        if(method=="rank-shared"||method=="rank-static")prefix_run<RankEngine>(input,method,swaps,maximum);
        else {require(method=="bitwise"||method=="peigen","method");prefix_run<Engine>(input,method,swaps,maximum);}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
