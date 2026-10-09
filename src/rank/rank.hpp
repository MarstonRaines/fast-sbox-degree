#pragma once
#include "degree.hpp"
#include "small_rank.hpp"

namespace coam {
// Coefficient u is the m-bit vector of coordinate-ANF coefficients at x^u.
// The histogram includes nonzero lambda only and uses deg(0)=0.
struct RankEngine {
    const Context& c;
    Box box, coefficients;
    Hist hist{}, old_hist{};
    uint64_t value=0, old_value=0;
    bool recompute, tiny;
    Word tiny_coefficients=0;
#ifdef COAM_COUNT
    uint64_t scanned_rows=0, basis_xors=0, coefficient_updates=0;
#endif

    RankEngine(const Context& context, const Box& input, bool from_scratch=false)
        : c(context), box(input), coefficients(input), recompute(from_scratch), tiny(c.n<=4 && c.m<=4) {
        if (box.size()!=c.length || std::any_of(box.begin(),box.end(),
                [&](uint32_t x){return x>=c.components;}))
            throw std::runtime_error("invalid rank-evaluator LUT");
        if (tiny) {
            tiny_coefficients=small_rank::transform(small_rank::pack(box.data(),c.n),c.n);
            Box().swap(coefficients);
        } else transform();
        profile();
    }

    // Materialize coefficients only for independent verification and inspection.
    Box coefficient_snapshot() const {
        if (!tiny) return coefficients;
        Box a(c.length);
        for (unsigned u=0;u<c.length;++u) a[u]=(tiny_coefficients>>(4*u))&15;
        return a;
    }

    void transform() {
        for (unsigned step=1;step<c.length;step*=2)
            for (unsigned base=0;base<c.length;base+=2*step)
                for (unsigned j=0;j<step;++j)
                    coefficients[base+step+j]^=coefficients[base+j];
    }

    void profile() {
        if (tiny) {value=small_rank::profile(tiny_coefficients,c.n,c.m,hist);return;}
        hist.fill(0);
        std::array<uint32_t,32> basis{};
        unsigned rank=0;
        uint32_t previous=c.components-1;
        size_t pos=0;
        for (unsigned d=c.n;d>0;--d) {
            while (pos<c.order.size() && c.weights[c.order[pos]]==d) {
                uint32_t v=coefficients[c.order[pos++]];
#ifdef COAM_COUNT
                ++scanned_rows;
#endif
                while (v) {
                    unsigned pivot=std::bit_width(v)-1;
                    if (!basis[pivot]) {basis[pivot]=v;++rank;break;}
                    v^=basis[pivot];
#ifdef COAM_COUNT
                    ++basis_xors;
#endif
                }
                if (rank==c.m) break;
            }
            uint32_t cumulative=(uint32_t(1)<<(c.m-rank))-1;
            hist[d]=previous-cumulative;
            previous=cumulative;
            if (rank==c.m) break; // All remaining lower-degree counts are zero.
        }
        hist[0]=previous;
        value=spectrum_sum(hist);
    }

    // Enumerate exactly A(p) symmetric-difference A(q), in 64-index blocks.
    // Block skipping uses the same supersets as the bit-packed method.
    void change_coefficients(unsigned p,unsigned q,uint32_t output_difference) {
        if (!output_difference || p==q) return;
        if (tiny) {
            tiny_coefficients^=(small_rank::tables.supersets[p]^small_rank::tables.supersets[q])*output_difference;
            return;
        }
        auto block=[&](unsigned w) {
            Word a=((w&(p>>6))==(p>>6))?c.supermask[p&63]:0;
            Word b=((w&(q>>6))==(q>>6))?c.supermask[q&63]:0;
            Word bits=a^b;
            if (c.n<6) bits&=(Word(1)<<c.length)-1;
            while (bits) {
                unsigned u=64*w+std::countr_zero(bits);
                coefficients[u]^=output_difference;
                bits&=bits-1;
#ifdef COAM_COUNT
                ++coefficient_updates;
#endif
            }
        };
        unsigned ph=p>>6,qh=q>>6;
        if (c.words>=16 && (c.words>>std::popcount(ph))+
                          (c.words>>std::popcount(qh))<c.words) {
            for (unsigned w=ph;w<c.words;w=(w+1)|ph) block(w);
            for (unsigned w=qh;w<c.words;w=(w+1)|qh)
                if ((w&ph)!=ph) block(w);
        } else for (unsigned w=0;w<c.words;++w) block(w);
    }

    uint64_t apply(unsigned p,unsigned q,bool reversible) {
        if (p>=c.length || q>=c.length) throw std::runtime_error("invalid swap");
        if (reversible) {old_hist=hist;old_value=value;}
        uint32_t difference=box[p]^box[q];
        std::swap(box[p],box[q]);
        if (!difference) return value;
        if (recompute) {
            if (tiny) tiny_coefficients=small_rank::transform(small_rank::pack(box.data(),c.n),c.n);
            else {coefficients=box;transform();}
        }
        else change_coefficients(p,q,difference);
        profile();return value;
    }

    // Supports immediate rollback of the most recent reversible trial.
    void undo(unsigned p,unsigned q) {
        if (recompute) {
            std::swap(box[p],box[q]);
            if (tiny) tiny_coefficients=small_rank::transform(small_rank::pack(box.data(),c.n),c.n);
            else {coefficients=box;transform();}
        }
        else {change_coefficients(p,q,box[p]^box[q]);std::swap(box[p],box[q]);}
        hist=old_hist;value=old_value;
    }
};
}
