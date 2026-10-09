#pragma once
#include <array>
#include <bit>
#include <cstdint>

namespace coam::small_min5 {
// Every linear subspace of F_2^5. Transitions depend only on a coefficient
// vector and its accumulated span, never on a benchmark instance.
struct Tables {
    std::array<uint32_t,374> spans{};
    std::array<std::array<uint16_t,32>,374> next{};
    std::array<uint8_t,374> ranks{};
    unsigned count=1;
    constexpr Tables() {
        std::array<unsigned,1024> index{};
        auto hash=[](uint32_t span){return (span*2654435761U)>>22;};
        spans[0]=1; index[hash(1)]=1;
        for(unsigned s=0;s<count;++s) {
            ranks[s]=uint8_t(std::bit_width(unsigned(std::popcount(spans[s])))-1);
            for(unsigned v=0;v<32;++v) {
                uint32_t enlarged=spans[s];
                for(unsigned x=0;x<32;++x)
                    if((spans[s]>>x)&1U)enlarged|=1U<<(x^v);
                unsigned slot=hash(enlarged);
                while(index[slot] && spans[index[slot]-1]!=enlarged)slot=(slot+1)&1023;
                if(!index[slot]){spans[count]=enlarged;index[slot]=++count;}
                next[s][v]=uint16_t(index[slot]-1);
            }
        }
    }
};
inline constexpr Tables tables{};
static_assert(tables.count==374);
using Packed=std::array<uint64_t,4>;
inline Packed pack(const uint32_t* box) {
    Packed a{};
    for(unsigned u=0;u<32;++u)a[u/8]|=uint64_t(box[u])<<(8*(u%8));
    return a;
}
inline Packed transform(Packed a) {
    for(auto&v:a){v^=(v<<8)&0xff00ff00ff00ff00ULL;v^=(v<<16)&0xffff0000ffff0000ULL;v^=v<<32;}
    a[1]^=a[0];a[3]^=a[2];a[2]^=a[0];a[3]^=a[1];
    return a;
}
inline constexpr std::array<unsigned,32> order{31,30,29,27,23,15,28,26,25,22,21,19,14,13,11,7,24,20,18,17,12,10,9,6,5,3,16,8,4,2,1,0};
inline unsigned minimum(const Packed& a,unsigned m) {
    unsigned state=0;
    for(unsigned u:order){
        state=tables.next[state][(a[u/8]>>(8*(u%8)))&31];
        if(tables.ranks[state]==m)return std::popcount(u);
    }
    return 0;
}
}
