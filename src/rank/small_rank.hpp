#pragma once
#include <array>
#include <bit>
#include <cstdint>

namespace coam::small_rank {
// All 67 linear subspaces of F_2^4, generated algebraically at compile time.
// A transition inserts one coefficient row. No input-specific data are stored.
struct Tables {
    std::array<uint16_t,67> spans{};
    std::array<std::array<uint8_t,16>,67> next{};
    std::array<uint8_t,67> ranks{};
    std::array<uint64_t,16> supersets{};
    unsigned count=1;
    constexpr Tables() {
        spans[0]=1;
        for (unsigned s=0;s<count;++s) {
            ranks[s]=uint8_t(std::bit_width(unsigned(std::popcount(spans[s])))-1);
            for (unsigned v=0;v<16;++v) {
                uint16_t enlarged=spans[s];
                for (unsigned x=0;x<16;++x)
                    if ((spans[s]>>x)&1U) enlarged|=uint16_t(1U<<(x^v));
                unsigned t=0;
                while (t<count && spans[t]!=enlarged) ++t;
                if (t==count) spans[count++]=enlarged;
                next[s][v]=uint8_t(t);
            }
        }
        for (unsigned p=0;p<16;++p)
            for (unsigned u=0;u<16;++u)
                if ((u&p)==p) supersets[p]|=uint64_t(1)<<(4*u);
    }
};
inline constexpr Tables tables{};
static_assert(tables.count==67);
inline uint64_t pack(const uint32_t* box,unsigned n) {
    uint64_t result=0;
    for (unsigned u=0;u<(1U<<n);++u) result|=uint64_t(box[u])<<(4*u);
    return result;
}
inline uint64_t transform(uint64_t a,unsigned n) {
    a^=(a<<4)&0xf0f0f0f0f0f0f0f0ULL;
    if(n>=2)a^=(a<<8)&0xff00ff00ff00ff00ULL;
    if(n>=3)a^=(a<<16)&0xffff0000ffff0000ULL;
    if(n>=4)a^=a<<32;
    return a;
}
inline constexpr std::array<unsigned,16> order4{15,14,13,11,7,12,10,9,6,5,3,8,4,2,1,0};
inline constexpr std::array<unsigned,8> order3{7,6,5,3,4,2,1,0};
inline unsigned minimum(uint64_t a,unsigned n,unsigned m) {
    unsigned state=0;
    const unsigned* order=n==4?order4.data():order3.data();
    if(n<3) {
        for(int d=int(n);d>=0;--d)for(unsigned u=0;u<(1U<<n);++u) {
            if(std::popcount(u)!=d)continue;
            state=tables.next[state][(a>>(4*u))&15];
            if(tables.ranks[state]==m)return unsigned(d);
        }
    } else for(unsigned j=0;j<(1U<<n);++j) {
        unsigned u=order[j];state=tables.next[state][(a>>(4*u))&15];
        if(tables.ranks[state]==m)return std::popcount(u);
    }
    return 0;
}
template<class Histogram> inline uint64_t profile(uint64_t a,unsigned n,unsigned m,Histogram& hist) {
    hist.fill(0);unsigned state=0,previous=(1U<<m)-1;uint64_t value=0;
    const unsigned* order=n==4?order4.data():order3.data();
    unsigned pos=0;
    for(unsigned d=n;d>0;--d) {
        if(n<3) {
            for(unsigned u=0;u<(1U<<n);++u)if(std::popcount(u)==int(d))
                state=tables.next[state][(a>>(4*u))&15];
        } else while(pos<(1U<<n) && std::popcount(order[pos])==int(d)) {
            unsigned u=order[pos++];state=tables.next[state][(a>>(4*u))&15];
            if(tables.ranks[state]==m)break;
        }
        unsigned cumulative=(1U<<(m-tables.ranks[state]))-1;
        hist[d]=previous-cumulative;value+=uint64_t(d)*hist[d];previous=cumulative;
        if(tables.ranks[state]==m)break;
    }
    hist[0]=previous;return value;
}
}
