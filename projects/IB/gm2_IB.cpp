#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <common.hpp>

size_t T=96;
size_t nm=3,nr=2;
int ib=2;
double a;
double mass_muon=0.1056583745;
index_t<4> ind;

//! read a single vector, for a specific mass and r, real or imaginary
djvec_t read(const char *what,size_t im,size_t r,size_t reim)
{return read_djvec(combine("data/corr%s",what),T,ind({im,im,r,reim}));}

//! read a combination of r and return appropriately simmetrized
djvec_t read(const char *what,size_t tpar,size_t im,int rpar,size_t reim)
{
  djvec_t o(T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,im,r,reim);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read VV averaging the three channels
djvec_t read_VV(const char *what,int tpar,size_t im,int rpar,size_t reim)
{
  return djvec_t(read(combine("%s_V1V1",what).c_str(),tpar,im,rpar,reim)+
		 read(combine("%s_V2V2",what).c_str(),tpar,im,rpar,reim)+
		 read(combine("%s_V3V3",what).c_str(),tpar,im,rpar,reim))/3.0;
}

//! compute the kernel f(Q2)
double kern_Q2(double Q2)
{
  double am=mass_muon*a;
  cout<<"a "<<a<<" "<<am<<endl;
  double s=Q2/sqr(am);
  //cout<<"s "<<s<<endl;
  double Z=(sqrt(1+4/s)-1)/2;
  //cout<<"Z "<<Z<<endl;
  return s*Z*Z*Z*(1-s*Z)/(1+s*Z*Z)/sqr(am);
}

//! compute the kernel f(t)
double kern_t(size_t t)
{
  double out=0;
  for(size_t iq=1;iq<T;iq++)
    {
      double Q=iq*2*M_PI/T;
      //cout<<"Q "<<Q<<endl;
      double Q2=sqr(Q);
      out+=2*kern_Q2(Q2)*((cos(Q*t)-1)/Q2+t*t/2);
  }
  return out;
}

int main()
{
  init_common_IB("/home/francesco/QCD/LAVORI/IB_MES_NF211/one_source/ultimate_input.txt");
  a=(1/lat_par[0].ainv[0]).ave();
  
  for(size_t t=0;t<T;t++) cout<<t<<" "<<kern_t(t)<<endl;
  return 0;
  ind.set_ranges({nm,nm,nr,2});
  
  const int il=0;
  djvec_t V0P5_LL=read("LL_V0P5",-1,il,-1,IM);
  djvec_t V0P5_0M=read("0M_V0P5",-1,il,-1,IM);
  djvec_t V0P5_0T=read("0T_V0P5",-1,il,-1,IM);
  djvec_t V0P5_0P=read("0P_V0P5",-1,il,-1,RE);
  
  djvec_t num_deltam_cr=forward_derivative(djvec_t(V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write("plots/num_deltam_cr.xmg");
  djvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write("plots/den_deltam_cr.xmg");
  
  djvec_t deltam_cr=-num_deltam_cr/den_deltam_cr;
  deltam_cr.ave_err().write("plots/deltam_cr_t.xmg");
  
  const int im=1;
  djvec_t VV_00=read_VV("00",1,im,1,RE);
  djvec_t VV_0P=read_VV("0P",1,im,-1,IM);
  djvec_t VV_0T=read_VV("0T",1,im,1,RE);
  djvec_t VV_0M=read_VV("0M",1,im,1,RE);
  djvec_t VV_LL=read_VV("LL",1,im,1,RE);
  effective_mass(VV_00).ave_err().write("plots/VV_00.xmg");
  djvec_t(VV_0P/VV_00).ave_err().write("plots/VV_0P.xmg");
  djvec_t(VV_0T/VV_00).ave_err().write("plots/VV_0T.xmg");
  djvec_t(VV_0M/VV_00).ave_err().write("plots/VV_0M.xmg");
  djvec_t(VV_LL/VV_00).ave_err().write("plots/VV_LL.xmg");
  
  djvec_t c=2.0*djvec_t(VV_0T+VV_0M)+VV_LL;
  c=c.subset(0,T/2-1);
  djvec_t d=deltam_cr*VV_0P.subset(0,T/2-1);
  c.ave_err().write("plots/c.xmg");
  d.ave_err().write("plots/d.xmg");
  
  return 0;
}
