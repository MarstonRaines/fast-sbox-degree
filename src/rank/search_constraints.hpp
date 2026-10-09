#pragma once
#include "common.hpp"

// Walsh/DDT transposition updates: A. V. Menyachikhin, Mathematical Aspects
// of Cryptography 11(2), 2020, 111--123, DOI 10.4213/mvk325.
// Walsh coefficients are unnormalised signed sums. All search arms share this.
struct SearchConstraints {
    unsigned n,size;
    std::vector<int16_t> walsh,ddt;
    std::vector<unsigned> walsh_counts,ddt_counts;
    std::vector<uint8_t> parity;
    unsigned walsh_max=0,ddt_max=0;
    SearchConstraints(unsigned bits,const Box& box)
        : n(bits),size(1U<<bits),walsh(size*size),ddt(size*size),
          walsh_counts(size+1),ddt_counts(size+1),parity(size) {
        require(n>=2&&n<=8&&box.size()==size,"constraint dimensions");
        for(unsigned x=0;x<size;++x) {
            require(box[x]<size,"constraint LUT entry");parity[x]=std::popcount(x)&1U;
        }
        for(unsigned l=1;l<size;++l) {
            auto* row=walsh.data()+l*size;
            for(unsigned x=0;x<size;++x)row[x]=1-2*int(parity[l&box[x]]);
            for(unsigned step=1;step<size;step*=2)
                for(unsigned base=0;base<size;base+=2*step)
                    for(unsigned j=0;j<step;++j) {
                        int a=row[base+j],b=row[base+step+j];
                        row[base+j]=a+b;row[base+step+j]=a-b;
                    }
            for(unsigned a=0;a<size;++a)++walsh_counts[std::abs(row[a])];
        }
        for(unsigned a=1;a<size;++a) {
            auto* row=ddt.data()+a*size;
            for(unsigned x=0;x<size;++x)++row[box[x]^box[x^a]];
            for(unsigned b=0;b<size;++b)++ddt_counts[row[b]];
        }
        refresh_max(walsh_counts,walsh_max);refresh_max(ddt_counts,ddt_max);
    }
    static void refresh_max(const std::vector<unsigned>& counts,unsigned& maximum) {
        maximum=counts.size()-1;while(maximum&&!counts[maximum])--maximum;
    }
    int nonlinearity() const {return int(size/2)-int(walsh_max/2);}
    unsigned uniformity() const {return ddt_max;}
    // Algorithm 1 (p. 114), with W=2^n LAT: each changed entry moves by +/-4.
    // old_p/old_q are the pre-update outputs. Reversing them undoes the update.
    void menyachikhin_walsh_update(unsigned p,unsigned q,unsigned old_p,unsigned old_q) {
        const unsigned delta=old_p^old_q,input_delta=p^q;
        if(!delta||!input_delta)return;
        std::array<unsigned,128> masks{};std::array<int,128> changes{};unsigned used=0;
        for(unsigned a=0;a<size;++a)if(parity[a&input_delta]) {
            masks[used]=a;changes[used++]=parity[a&p]?4:-4;
        }
        for(unsigned l=1;l<size;++l)if(parity[l&delta]) {
            auto* row=walsh.data()+l*size;int sign=parity[l&old_p]?-1:1;
            for(unsigned j=0;j<used;++j) {
                auto& v=row[masks[j]];unsigned before=std::abs(v);
                v+=sign*changes[j];unsigned after=std::abs(v);
                --walsh_counts[before];++walsh_counts[after];
            }
        }
        refresh_max(walsh_counts,walsh_max);
    }
    void shift_ddt(unsigned a,unsigned b,int amount) {
        auto& v=ddt[a*size+b];--ddt_counts[v];v+=amount;++ddt_counts[v];
    }
    // Algorithm 2 (p. 117), in integer DDT counts (2^n times probabilities).
    // LUT is AFTER swapping p,q. a=p xor q is unchanged; the other rows have
    // two changed unordered input pairs, each contributing twice. Undo while
    // the candidate LUT is still present, then restore the LUT in the caller.
    void menyachikhin_ddt_update(const Box& after,unsigned p,unsigned q,bool undo=false) {
        if(p==q||after[p]==after[q])return;
        int amount=undo?-2:2;
        for(unsigned a=1;a<size;++a)if(a!=(p^q)) {
            unsigned xp=after[p^a],xq=after[q^a];
            unsigned b0=after[q]^xp,b2=after[p]^xq,b1=after[p]^xp;
            if(b0==b2) { // Algorithm 2's coincident-entry case, +/-4.
                shift_ddt(a,b0,-2*amount);shift_ddt(a,b1,2*amount);
            } else { // Four entries, alternating -2,+2,-2,+2.
                unsigned b3=after[q]^xq;
                shift_ddt(a,b0,-amount);shift_ddt(a,b1,amount);
                shift_ddt(a,b2,-amount);shift_ddt(a,b3,amount);
            }
        }
        refresh_max(ddt_counts,ddt_max);
    }
    void verify(const Box& box,bool check_walsh=true,bool check_ddt=true) const {
        SearchConstraints reference(n,box);
        if(check_walsh)require(walsh==reference.walsh&&walsh_counts==reference.walsh_counts&&
            walsh_max==reference.walsh_max,"incremental Walsh state/max mismatch");
        if(check_ddt)require(ddt==reference.ddt&&ddt_counts==reference.ddt_counts&&
            ddt_max==reference.ddt_max,"incremental DDT state/max mismatch");
    }
};
