#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <common.hpp>
#include <gm2_IB_integrators.hpp>

const int ilight=0;
size_t T=48;//96;
size_t nm=3,nr=2;
int ib=0;//2;
const int ens_id=3;
index_t<4> ind;
double Za=0.731;//0.762;
boot_init_t bi;

const int tmin_cr=12,tmax_cr=22;
const int im=2;
const double eq=ec;

//! read a single vector, for a specific mass and r, real or imaginary
dbvec_t read(const char *what,size_t im,size_t r,size_t reim)
{return dbvec_t(bi,read_djvec(combine("data/corr%s",what),T,ind({im,im,r,reim})));}

//! read a combination of r and return appropriately simmetrized
dbvec_t read(const char *what,int tpar,size_t im,int rpar,size_t reim)
{
  dbvec_t o(T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,im,r,reim)*(r==0?1:rpar);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read VV averaging the three channels
dbvec_t read_VV_ren(const char *what,int tpar,size_t im,int rpar,size_t reim)
{
  return Za*Za*dbvec_t(read(combine("%s_V1V1",what).c_str(),tpar,im,rpar,reim)+
		       read(combine("%s_V2V2",what).c_str(),tpar,im,rpar,reim)+
		       read(combine("%s_V3V3",what).c_str(),tpar,im,rpar,reim))/3.0;
}

//! integrate
dboot_t integrate(const dbvec_t &corr,const dboot_t &a)
{
  dboot_t out;
  out=0.0;
  for(size_t t=1;t<T/2-4;t++) out+=corr[t]*ftilde_t(t,a);
  out*=4*sqr(alpha_em)*sqr(eq);
  
  return out;
}

//! compute the critical deltam
dboot_t compute_deltam_cr(size_t iq)
{
  dbvec_t V0P5_LL=read("LL_V0P5",-1,iq,-1,IM);
  dbvec_t V0P5_0M=read("0M_V0P5",-1,iq,-1,IM);
  dbvec_t V0P5_0T=read("0T_V0P5",-1,iq,-1,IM);
  dbvec_t num_deltam_cr=forward_derivative(dbvec_t(V0P5_LL+2.0*dbvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write("plots/num_deltam_cr.xmg");

  dbvec_t V0P5_0P=read("0P_V0P5",-1,iq,+1,RE);
  dbvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write("plots/den_deltam_cr.xmg");
  
  dboot_t deltam_cr=constant_fit(dbvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),tmin_cr,tmax_cr,"plots/deltam_cr_t.xmg");
  
  return deltam_cr;
}

//! load VV for the leading order
dbvec_t load_LO(const int im)
{return read_VV_ren("00",1,im,1,RE);}

//! load QED corrections
dbvec_t load_QED(const int im,const dboot_t &deltam_cr,const dbvec_t &VV_LO)
{
  dbvec_t VV_0T=read_VV_ren("0T",1,im,1,RE);
  dbvec_t VV_0M=read_VV_ren("0M",1,im,1,RE);
  dbvec_t VV_LL=read_VV_ren("LL",1,im,1,RE);
  dbvec_t(VV_0T/VV_LO).ave_err().write("plots/VV_0T.xmg");
  dbvec_t(VV_0M/VV_LO).ave_err().write("plots/VV_0M.xmg");
  dbvec_t(VV_LL/VV_LO).ave_err().write("plots/VV_LL.xmg");
  dbvec_t c=dbvec_t(VV_LL+2.0*dbvec_t(VV_0T+VV_0M));
  
  dbvec_t VV_0P=read_VV_ren("0P",1,im,-1,IM);
  dbvec_t(VV_0P/VV_LO).ave_err().write("plots/VV_0P.xmg");
  dbvec_t d=-(deltam_cr*VV_0P);
  return dbvec_t(c+2.0*d)*e2*sqr(eq);
}

int main()
{
  init_common_IB("/home/francesco/QCD/LAVORI/IB_MES_NF211/one_source/ultimate_input.txt");
  dboot_t a=1/lat_par[0].ainv[ib];
  int input_an_id=0;
  bi=jack_index[input_an_id][ens_id];
  cout.precision(16);
  
  // for(size_t t=1;t<1000;t++) cout<<t<<" "<<(ftilde_t(t,a)/sqr(t)).ave_err()<<endl;
  // return 0;
  ind.set_ranges({nm,nm,nr,2});
  
  dboot_t deltam_cr=compute_deltam_cr(ilight);
  
  dbvec_t VV_LO=load_LO(im);
  dbvec_t VV_QED=load_QED(im,deltam_cr,VV_LO);
  //effective_mass(VV_LO).ave_err().write("plots/VV_LO.xmg");
  
  // dbvec_t c=2.0*dbvec_t(VV_0T+VV_0M)+VV_LL;
  // c=c.subset(0,T/2-1);
  // dbvec_t d=deltam_cr*VV_0P.subset(0,T/2-1);
  // c.ave_err().write("plots/c.xmg");
  // d.ave_err().write("plots/d.xmg");
  
  dboot_t amu=integrate(VV_LO,a);
  cout<<"amu: "<<amu.ave_err()<<endl;
  
  dboot_t amu_QED=integrate(VV_QED,a);
  cout<<"amu_QED: "<<amu_QED.ave_err()<<endl;
  
  cout<<" Ratio: "<<dboot_t(amu_QED/amu).ave_err()<<endl;
  
  close_integrators();
  
  return 0;
}
