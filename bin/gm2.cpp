#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

size_t T=48;
size_t nm=1,nr=1;
index_t<4> ind;

djvec_t read(const char *what,size_t tpar,size_t reim)
{return read_djvec(combine("data/mes_contr_%s",what),T,ind({0,0,0,reim})).symmetrized(tpar);}

djvec_t read_VV(const char *what,size_t tpar,size_t reim)
{
  return djvec_t(read(combine("V1V1_%s",what).c_str(),tpar,reim)+
		 read(combine("V2V2_%s",what).c_str(),tpar,reim)+
		 read(combine("V3V3_%s",what).c_str(),tpar,reim))/3.0;
}

int main()
{
  set_njacks(16);
  ind.set_ranges({nm,nm,nr,2});
  
  djvec_t V0P5_LL=read("V0P5_LL",-1,IM);
  djvec_t V0P5_0M=read("V0P5_0M",-1,IM);
  djvec_t V0P5_0T=read("V0P5_0T",-1,IM);
  djvec_t V0P5_0P=read("V0P5_0P",-1,RE);
  
  djvec_t num_deltam_cr=forward_derivative(djvec_t(V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write("plots/num_deltam_cr.xmg");
  djvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write("plots/den_deltam_cr.xmg");
  
  djvec_t deltam_cr=-num_deltam_cr/den_deltam_cr;
  deltam_cr.ave_err().write("plots/deltam_cr_t.xmg");
  
  djvec_t VV_00=read_VV("00",1,RE);
  djvec_t VV_0P=read_VV("0P",1,IM);
  djvec_t VV_0T=read_VV("0T",1,RE);
  djvec_t VV_0M=read_VV("0M",1,RE);
  djvec_t VV_LL=read_VV("LL",1,RE);
  VV_00.ave_err().write("plots/VV_00.xmg");
  VV_0P.ave_err().write("plots/VV_0P.xmg");
  VV_0T.ave_err().write("plots/VV_0T.xmg");
  VV_0M.ave_err().write("plots/VV_0M.xmg");
  VV_LL.ave_err().write("plots/VV_LL.xmg");
  
  djvec_t c=2.0*djvec_t(VV_0T+VV_0M)+VV_LL;
  c=c.subset(0,T/2-1);
  djvec_t d=deltam_cr*VV_0P.subset(0,T/2-1);
  c.ave_err().write("plots/c.xmg");
  d.ave_err().write("plots/d.xmg");
  
  return 0;
}
