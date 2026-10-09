#include "common.hpp"
void verify(const Context& c,const Box& box,unsigned& states) {
    Engine mini(c,Method::Proposed,Metric::Minimum,box,true),maximum(c,Method::Proposed,Metric::Maximum,box,true);
    auto check=[&] {
        auto h=independent_histogram(c,mini.box);
        require(mini.value==hist_min(h)&&maximum.value==hist_max(h),"minimum/maximum mismatch");
        auto coefficients=independent_coefficients(mini.box);std::vector<unsigned> expected;
        for(unsigned l=1;l<c.components;++l){unsigned degree=0;
            for(unsigned u=0;u<c.length;++u)if(std::popcount(l&coefficients[u])&1U)degree=std::max(degree,unsigned(std::popcount(u)));
            if(degree==mini.value)expected.push_back(l);
        }
        require(mini.minimum_lambdas()==expected,"minimum-vector enumeration mismatch");++states;
    };
    check();
    for(auto[p,q]:std::array<Swap,3>{{{0,1},{0,c.length-1},{1,c.length-1}}}) {
        auto prior=mini.box;mini.apply(p,q,true);maximum.apply(p,q,true);check();
        mini.undo(p,q);maximum.undo(p,q);mini.value=mini.evaluate();
        require(mini.box==prior&&maximum.box==prior,"rollback input mismatch");check();
    }
}
void verify_minimum_updates(const Context& c,const Box& box,unsigned& states) {
    for(bool prepared:{false,true}) {
        Engine mini(c,Method::Proposed,Metric::Minimum,box,prepared);
        auto check=[&] {
            auto expected=hist_min(independent_histogram(c,mini.box));
            require(mini.value==expected,"updated minimum mismatch");
            require(mini.evaluate()==expected,"cached minimum mismatch");
            Engine fresh(c,Method::Proposed,Metric::Minimum,mini.box,prepared);
            require(fresh.value==expected,"fresh minimum mismatch");++states;
        };
        constexpr std::array<Swap,5> swaps{{{0,1},{4,6},{8,12},{16,24},{2,18}}};
        check();
        for(auto[p,q]:swaps) {
            auto prior=mini.box;mini.apply(p,q,true);check();
            mini.undo(p,q);require(mini.box==prior,"minimum rollback input mismatch");check();
            mini.apply(p,q,false);check();
        }
        for(auto[p,q]:swaps) {
            auto prior=mini.box;mini.apply(p,q,true);check();
            mini.undo(p,q);require(mini.box==prior,"minimum rollback input mismatch");check();
        }
        auto prior=mini.box;mini.apply(0,0,true);check();
        mini.undo(0,0);require(mini.box==prior,"minimum self-swap input mismatch");check();
    }
}
int main(){try{unsigned states=0;std::mt19937 random(20261005);
    for(unsigned n=1;n<=7;++n)for(unsigned m=1;m<=n;++m){Context c(n,m);
        for(unsigned kind=0;kind<6;++kind){Box box(c.length);
            for(unsigned x=0;x<c.length;++x)box[x]=kind==0?0:kind==1?c.components-1:kind==2?x&(c.components-1):random()%c.components;
            verify(c,box,states);
            if(n==5)verify_minimum_updates(c,box,states);
        }}
    std::cout<<"{\"passed\":true,\"states\":"<<states<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
