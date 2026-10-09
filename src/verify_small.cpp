#include "rank/small_rank.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>
using namespace coam;
using Box=std::vector<uint32_t>;
using Hist=std::array<uint32_t,21>;
uint64_t states=0,swaps_checked=0;
void require(bool p){if(!p)throw std::runtime_error("small-path oracle mismatch");}
Hist oracle(const Box& s,unsigned n,unsigned m){
    Hist h{};
    for(unsigned lambda=1;lambda<(1U<<m);++lambda){
        unsigned degree=0;
        for(unsigned u=0;u<(1U<<n);++u){
            unsigned v=0,x=u;
            do {v^=s[x];if(!x)break;x=(x-1)&u;}while(true);
            if(std::popcount(v&lambda)%2)degree=std::max(degree,unsigned(std::popcount(u)));
        }
        ++h[degree];
    }
    return h;
}
void check(const Box& s,unsigned n,unsigned m,uint64_t a){
    Hist expected=oracle(s,n,m),actual{};uint64_t sum=0;
    for(unsigned d=0;d<=n;++d)sum+=d*expected[d];
    require(small_rank::profile(a,n,m,actual)==sum && actual==expected);
    unsigned minimum=0;while(!expected[minimum])++minimum;
    require(small_rank::minimum(a,n,m)==minimum);
    ++states;
}
void exercise(Box s,unsigned n,unsigned m,bool updates){
    auto a=small_rank::transform(small_rank::pack(s.data(),n),n);check(s,n,m,a);
    if(!updates)return;
    for(unsigned p=0;p<s.size();++p)for(unsigned q=p+1;q<s.size();++q){
        auto delta=(small_rank::tables.supersets[p]^small_rank::tables.supersets[q])*(s[p]^s[q]);
        std::swap(s[p],s[q]);check(s,n,m,a^delta);
        std::swap(s[p],s[q]);require(((a^delta)^delta)==a);++swaps_checked;
    }
}
int main(){
    for(unsigned s=0;s<67;++s)for(unsigned v=0;v<16;++v){
        auto span=small_rank::tables.spans[s];unsigned expected=span;
        for(unsigned x=0;x<16;++x)if((span>>x)&1U)expected|=1U<<(x^v);
        auto t=small_rank::tables.next[s][v];require(expected==small_rank::tables.spans[t]);
        require((1U<<small_rank::tables.ranks[t])==unsigned(std::popcount(expected)));
    }
    Box b(8);std::iota(b.begin(),b.end(),0);do{exercise(b,3,3,true);}while(std::next_permutation(b.begin(),b.end()));
    for(unsigned k=0;k<65536;++k){Box s(4);for(unsigned u=0;u<4;++u)s[u]=(k>>(4*u))&15;exercise(s,2,4,false);}
    std::mt19937 rng(20261007);
    for(unsigned n=1;n<=4;++n)for(unsigned m=1;m<=4;++m)for(unsigned k=0;k<256;++k){
        Box s(1U<<n);for(auto&v:s)v=rng()&((1U<<m)-1);exercise(s,n,m,true);
    }
    std::cout<<"{\"passed\":true,\"transition_checks\":1072,\"states\":"<<states<<",\"swaps\":"<<swaps_checked<<"}\n";
}
