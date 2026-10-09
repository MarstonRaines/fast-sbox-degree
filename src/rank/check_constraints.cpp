#include "search_constraints.hpp"
#include <numeric>

void direct_walsh_samples(const SearchConstraints& s,const Box& box) {
    unsigned stride=s.size<=16?1:17;
    for(unsigned l=1;l<s.size;l+=stride)for(unsigned a=0;a<s.size;a+=stride) {
        int value=0;
        for(unsigned x=0;x<s.size;++x)
            value+=1-2*int((std::popcount(l&box[x])+std::popcount(a&x))&1U);
        require(value==s.walsh[l*s.size+a],"direct Walsh sum mismatch");
    }
}
int main(int argc,char** argv) {
 try {
    { // Published Examples 1 and 2, including both DDT branches.
        Box g{2,4,1,0,3,5,6,7};SearchConstraints c(3,g);const auto old=c.walsh;
        c.menyachikhin_walsh_update(1,2,g[1],g[2]);std::swap(g[1],g[2]);c.verify(g,true,false);
        require(c.walsh[1*8+1]-old[1*8+1]==4&&c.walsh[3*8+2]-old[3*8+2]==-4,"Menyachikhin Example 1");
        Box h{5,6,1,7,2,4,3,0};SearchConstraints d(3,h);const auto old_ddt=d.ddt;
        std::swap(h[0],h[1]);d.menyachikhin_ddt_update(h,0,1);d.verify(h,false,true);
        require(d.ddt[4*8+4]-old_ddt[4*8+4]==2&&d.ddt[7*8+5]-old_ddt[7*8+5]==-4,"Menyachikhin Example 2");
    }
    std::mt19937_64 rng(0x20261008ULL);uint64_t swaps=0,commits=0,walsh_only=0;
    unsigned nl_up=0,nl_down=0,du_up=0,du_down=0;
    auto exercise=[&](Box& box,SearchConstraints& s,unsigned p,unsigned q,unsigned mode) {
        const auto before=s;const auto old=box;
        int old_nl=s.nonlinearity();unsigned old_du=s.uniformity();
        std::swap(box[p],box[q]);s.menyachikhin_walsh_update(p,q,old[p],old[q]);
        s.verify(box,true,false);direct_walsh_samples(s,box);
        nl_up+=s.nonlinearity()>old_nl;nl_down+=s.nonlinearity()<old_nl;
        if(mode==0) { // Rejection at the NL stage leaves DDT untouched.
            s.menyachikhin_walsh_update(p,q,old[q],old[p]);box=old;++walsh_only;
        } else {
            s.menyachikhin_ddt_update(box,p,q);s.verify(box);
            du_up+=s.uniformity()>old_du;du_down+=s.uniformity()<old_du;
            if(mode==1) {s.menyachikhin_ddt_update(box,p,q,true);s.menyachikhin_walsh_update(p,q,old[q],old[p]);box=old;}
            else ++commits;
        }
        s.verify(box);direct_walsh_samples(s,box);
        if(mode!=2)require(s.walsh==before.walsh&&s.ddt==before.ddt&&
            s.walsh_counts==before.walsh_counts&&s.ddt_counts==before.ddt_counts,"rollback snapshot mismatch");
        ++swaps;
    };
    for(unsigned n=2;n<=5;++n)for(unsigned shape=0;shape<3;++shape) {
        Box box(1U<<n);std::iota(box.begin(),box.end(),0);
        if(shape==1)std::shuffle(box.begin(),box.end(),rng);
        if(shape==2)for(auto& x:box)x=rng()%box.size();
        SearchConstraints s(n,box);s.verify(box);
        for(unsigned p=0;p<box.size();++p)for(unsigned q=p;q<box.size();++q)
            exercise(box,s,p,q,(p+q)%3);
    }
    Box random(256);std::iota(random.begin(),random.end(),0);std::shuffle(random.begin(),random.end(),rng);
    std::vector<Box> large{random};for(int i=1;i<argc;++i)large.push_back(read_input(argv[i]).box);
    for(auto box:large) {
        SearchConstraints s(8,box);
        for(unsigned i=0;i<256;++i)exercise(box,s,rng()%256,rng()%256,i%3);
    }
    require(nl_up&&nl_down&&du_up&&du_down,"both directions of extrema must be covered");
    std::cout<<"{\"passed\":true,\"swaps\":"<<swaps<<",\"commits\":"<<commits
      <<",\"walsh_only_rollbacks\":"<<walsh_only<<",\"nl_increases\":"<<nl_up<<",\"nl_decreases\":"<<nl_down
      <<",\"du_increases\":"<<du_up<<",\"du_decreases\":"<<du_down<<"}\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
