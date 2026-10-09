#include "rank/small_min5.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using namespace coam;
void require(bool p){if(!p)throw std::runtime_error("5-bit independent oracle mismatch");}
unsigned oracle(const std::vector<uint32_t>& box,unsigned m,const small_min5::Packed& a){
    unsigned result=5;
    std::array<unsigned,32> coefficients{};
    for(unsigned u=0;u<32;++u){unsigned x=u;do{coefficients[u]^=box[x];if(!x)break;x=(x-1)&u;}while(true);require(coefficients[u]==((a[u/8]>>(8*(u%8)))&255));}
    for(unsigned lambda=1;lambda<(1U<<m);++lambda){unsigned d=0;for(unsigned u=0;u<32;++u)if(std::popcount(lambda&coefficients[u])%2)d=std::max(d,unsigned(std::popcount(u)));result=std::min(d,result);}
    return result;
}
int main(){
    for(unsigned s=0;s<374;++s)for(unsigned v=0;v<32;++v){
        auto span=small_min5::tables.spans[s];auto expected=span;
        for(unsigned x=0;x<32;++x)if((span>>x)&1U)expected|=1U<<(x^v);
        auto t=small_min5::tables.next[s][v];require(expected==small_min5::tables.spans[t]);require((1U<<small_min5::tables.ranks[t])==unsigned(std::popcount(expected)));
    }
    std::mt19937 rng(20261007);uint64_t states=0,swaps=0;
    for(unsigned m=1;m<=5;++m)for(unsigned k=0;k<64;++k){
        std::vector<uint32_t> box(32);for(unsigned x=0;x<32;++x)box[x]=k==0?0:k==1?((1U<<m)-1):k==2?(x&((1U<<m)-1)):(rng()&((1U<<m)-1));
        auto packed=small_min5::pack(box.data());
        auto check=[&]{auto a=small_min5::transform(packed);require(small_min5::minimum(a,m)==oracle(box,m,a));++states;};check();
        for(unsigned p=0;p<32;++p)for(unsigned q=p+1;q<32;++q){auto before=packed;auto d=box[p]^box[q];packed[p/8]^=uint64_t(d)<<(8*(p%8));packed[q/8]^=uint64_t(d)<<(8*(q%8));std::swap(box[p],box[q]);check();std::swap(box[p],box[q]);packed[p/8]^=uint64_t(d)<<(8*(p%8));packed[q/8]^=uint64_t(d)<<(8*(q%8));require(packed==before);++swaps;}
    }
    std::cout<<"{\"passed\":true,\"transition_checks\":11968,\"states\":"<<states<<",\"swaps\":"<<swaps<<",\"table_bytes\":"<<sizeof(small_min5::Tables)<<"}\n";
}
