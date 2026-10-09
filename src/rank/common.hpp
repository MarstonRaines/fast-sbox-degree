#pragma once
#include "rank.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sys/resource.h>

using namespace coam;
using Clock=std::chrono::steady_clock;
using Swap=std::pair<unsigned,unsigned>;
inline void require(bool condition,const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
inline uint64_t elapsed(Clock::time_point a,Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b-a).count();
}
inline uint64_t fold(uint64_t h,uint64_t value) {return (h^value)*1099511628211ULL;}
inline uint64_t box_hash(const Box& box) {
    uint64_t result=14695981039346656037ULL;
    for (auto x:box) result=fold(result,x);
    return result;
}
template<class C> void json_array(const C& values) {
    std::cout<<'[';bool first=true;
    for(auto x:values){if(!first)std::cout<<',';first=false;std::cout<<+x;}
    std::cout<<']';
}
struct Input {unsigned n,m;Box box;};
inline Input read_input(const std::string& name) {
    std::ifstream file(name);Input input{};
    require(bool(file>>input.n>>input.m),"cannot read input: "+name);
    require(input.n>=1 && input.n<=20 && input.m>=1 && input.m<=20,"dimensions");
    input.box.resize(1U<<input.n);
    for(auto& v:input.box)require(bool(file>>v)&&v<(1U<<input.m),"invalid LUT entry");
    std::string extra;require(!(file>>extra),"extra LUT entry");return input;
}
inline std::vector<Swap> read_swaps(const std::string& name,unsigned length) {
    std::ifstream file(name);require(bool(file),"cannot read swaps");
    std::vector<Swap> result;unsigned p,q;
    while(file>>p){require(bool(file>>q)&&p<length&&q<length&&p!=q,"invalid swap");result.emplace_back(p,q);}
    require(file.eof()&&!result.empty(),"empty/invalid swaps");return result;
}
inline Hist independent_histogram(const Context& c,const Box& box) {
    Hist result{};std::vector<uint8_t> a(c.length);
    for(unsigned l=1;l<c.components;++l) {
        for(unsigned x=0;x<c.length;++x)a[x]=std::popcount(l&box[x])&1U;
        for(unsigned bit=1;bit<c.length;bit<<=1)
            for(unsigned u=0;u<c.length;++u)if(u&bit)a[u]^=a[u^bit];
        unsigned d=0;
        for(unsigned u=0;u<c.length;++u)if(a[u])d=std::max(d,unsigned(std::popcount(u)));
        ++result[d];
    }
    return result;
}
// Direct subset summation, independently of the transform implementation.
inline Box independent_coefficients(const Box& box) {
    Box a(box.size());
    for(unsigned u=0;u<box.size();++u) {
        unsigned x=u;
        for(;;){a[u]^=box[x];if(!x)break;x=(x-1)&u;}
    }
    return a;
}
inline void warm_cpu() {
    auto start=Clock::now();uint64_t x=1;
    do {for(unsigned i=0;i<10000;++i)x=fold(x,i);asm volatile("" : "+r"(x));}
    while(elapsed(start,Clock::now())<50000000ULL);
}
inline uint64_t peak_rss_bytes() {
#ifdef __linux__
    // ru_maxrss can retain the pre-exec Python launcher's high-water mark.
    // /proc/self/status measures this executable's address space after exec.
    std::ifstream status("/proc/self/status");std::string key,line;
    while(status>>key) {
        if(key=="VmHWM:"){uint64_t kib=0;status>>kib;return kib*1024;}
        std::getline(status,line);
    }
    throw std::runtime_error("cannot read VmHWM");
#else
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    return uint64_t(usage.ru_maxrss);
#else
    return uint64_t(usage.ru_maxrss)*1024;
#endif
#endif
}
template<class T> uint64_t storage(const std::vector<T>& v){return v.capacity()*sizeof(T);}
inline uint64_t context_storage(const Context& c) {
    return sizeof(c)+storage(c.order)+storage(c.weights)+storage(c.degree_masks);
}
inline uint64_t state_storage(const RankEngine& e) {
    return sizeof(e)+storage(e.box)+storage(e.coefficients);
}
inline uint64_t state_storage(const Engine& e) {
    return sizeof(e)+storage(e.box)+storage(e.coord)+storage(e.work)+
        storage(e.combination)+storage(e.truth)+storage(e.degrees)+storage(e.old_degrees)+
        storage(e.output_anfs)+storage(e.output_truth)+
        storage(e.delta);
}
