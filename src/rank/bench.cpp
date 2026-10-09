#include "common.hpp"

template<class E> void timed_updates(const std::string& method,const Input& input,
        const std::vector<Swap>& swaps,uint64_t steps) {
    warm_cpu();
    auto begin=Clock::now();Context c(input.n,input.m);auto context_end=Clock::now();
    std::unique_ptr<E> engine;
    if constexpr(std::is_same_v<E,RankEngine>)
        engine=std::make_unique<E>(c,input.box,method=="rank-static");
    else engine=std::make_unique<E>(c,method_of(method),Metric::Spectrum,input.box,true);
    auto init_end=Clock::now();uint64_t digest=14695981039346656037ULL;
    for(uint64_t i=0;i<steps;++i) {
        auto [p,q]=swaps[i%swaps.size()];
        digest=fold(digest,engine->apply(p,q,false));
    }
    auto end=Clock::now();const auto& e=*engine;
    std::cout<<"{\"method\":"<<std::quoted(method)<<",\"n\":"<<input.n<<",\"m\":"<<input.m
       <<",\"steps\":"<<steps<<",\"context_ns\":"<<elapsed(begin,context_end)
       <<",\"init_ns\":"<<elapsed(context_end,init_end)<<",\"kernel_ns\":"<<elapsed(init_end,end)
       <<",\"total_ns\":"<<elapsed(begin,end)<<",\"peak_rss_bytes\":"<<peak_rss_bytes()
       <<",\"context_allocated_bytes\":"<<context_storage(c)<<",\"state_allocated_bytes\":"<<state_storage(e)
       <<",\"digest\":"<<digest<<",\"final_box_hash\":"<<box_hash(e.box)
       <<",\"final_sum\":"<<e.value<<",\"histogram\":";json_array(e.hist);
#ifdef COAM_COUNT
    if constexpr(std::is_same_v<E,RankEngine>)
        std::cout<<",\"scanned_rows\":"<<e.scanned_rows<<",\"basis_xors\":"<<e.basis_xors
                 <<",\"coefficient_updates\":"<<e.coefficient_updates;
#endif
    std::cout<<"}\n";
}

void test_state(const Context& c,const Box& initial,const std::vector<Swap>& swaps,
                uint64_t& states,uint64_t& updates) {
    RankEngine e(c,initial),static_e(c,initial,true);
    auto check=[&] {
        auto oracle=independent_histogram(c,e.box);
        require(e.hist==oracle && static_e.hist==oracle,"independent histogram mismatch");
        require(e.coefficient_snapshot()==independent_coefficients(e.box),"independent coefficient mismatch");
        require(static_e.coefficient_snapshot()==e.coefficient_snapshot(),"static coefficient mismatch");
        uint64_t count=0;for(auto v:e.hist)count+=v;
        require(count==c.components-1,"histogram total");++states;
    };
    check();
    for(size_t i=0;i<swaps.size();++i) {
        auto [p,q]=swaps[i];auto old_box=e.box,old_a=e.coefficient_snapshot();auto old_h=e.hist;
        e.apply(p,q,true);static_e.apply(p,q,true);++updates;check();
        if(i%3!=0){e.undo(p,q);static_e.undo(p,q);require(e.box==old_box&&e.coefficient_snapshot()==old_a&&e.hist==old_h,"rollback");check();}
    }
}

void check_all() {
    uint64_t maps=0,states=0,updates=0;std::mt19937 random(20261005);
    Context tiny(2,2);std::vector<Swap> all;
    for(unsigned p=0;p<4;++p)for(unsigned q=p;q<4;++q)all.emplace_back(p,q);
    for(unsigned encoded=0;encoded<256;++encoded){Box b(4);for(unsigned x=0;x<4;++x)b[x]=(encoded>>(2*x))&3;test_state(tiny,b,all,states,updates);++maps;}
    for(unsigned n=1;n<=8;++n)for(unsigned m=1;m<=n;++m) {
        Context c(n,m);std::vector<Swap> pairs;
        for(unsigned j=0;j<16;++j)pairs.emplace_back(random()%c.length,random()%c.length);
        for(unsigned kind=0;kind<8;++kind) {
            Box b(c.length);
            for(unsigned x=0;x<c.length;++x) {
                if(kind==0)b[x]=0;
                else if(kind==1)b[x]=c.components-1;
                else if(kind==2)b[x]=x&(c.components-1);
                else if(kind==3)b[x]=((x&3)==3)?c.components-1:0;
                else b[x]=random()%c.components;
            }
            test_state(c,b,pairs,states,updates);++maps;
        }
    }
    std::cout<<"{\"status\":\"pass\",\"maps\":"<<maps<<",\"states\":"<<states
             <<",\"updates\":"<<updates<<",\"independent_oracles\":[\"component_byte_transform\",\"direct_subset_sum\"]}\n";
}

void cross_check(const Input& input,const std::vector<Swap>& swaps,unsigned steps) {
    Context c(input.n,input.m);RankEngine r(c,input.box);
    Engine bit(c,Method::Bitwise,Metric::Spectrum,input.box,true);
    auto check=[&]{require(r.hist==bit.hist,"cross histogram");
        std::vector<Word> coord(c.m*c.words);c.coordinates(r.box,coord);
        auto coefficients=r.coefficient_snapshot();
        for(unsigned u=0;u<c.length;++u)for(unsigned b=0;b<c.m;++b)
            require(((coefficients[u]>>b)&1U)==((coord[b*c.words+u/64]>>(u%64))&1U),"coordinate recomputation");};
    check();
    for(unsigned i=0;i<steps;++i){auto[p,q]=swaps[i%swaps.size()];r.apply(p,q,true);bit.apply(p,q,true);check();
        if(i%3!=0){r.undo(p,q);bit.undo(p,q);check();}}
    std::cout<<"{\"status\":\"pass\",\"n\":"<<input.n<<",\"m\":"<<input.m<<",\"updates\":"<<steps<<"}\n";
}

int main(int argc,char** argv) {
    try {
        require(argc>=2,"check | cross INPUT SWAPS STEPS | bench METHOD INPUT SWAPS STEPS");
        std::string mode=argv[1];
        if(mode=="check"){check_all();return 0;}
        if(mode=="cross"){require(argc==5,"cross INPUT SWAPS STEPS");auto input=read_input(argv[2]);cross_check(input,read_swaps(argv[3],input.box.size()),std::stoul(argv[4]));return 0;}
        require(mode=="bench"&&argc==6,"bench METHOD INPUT SWAPS STEPS");
        auto input=read_input(argv[3]);auto swaps=read_swaps(argv[4],input.box.size());uint64_t steps=std::stoull(argv[5]);
        std::string method=argv[2];
        if(method=="rank-shared"||method=="rank-static")timed_updates<RankEngine>(method,input,swaps,steps);
        else {require(method=="bitwise"||method=="peigen","method");timed_updates<Engine>(method,input,swaps,steps);}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
