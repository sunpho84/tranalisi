#ifndef _GM2_IB_COMMON_HPP
#define _GM2_IB_COMMON_HPP

#include <common.hpp>
#include <gm2_IB_integrators.hpp>

const size_t ilight=0,istrange=1,icharm=2;
const vector<string> qname({"light","strange","charm"});
const double M_V_phys[3]={0.775,1.0195,3.0969};
size_t nm,nr;
index_t<4> ind;
const vector<vector<ave_err_t>> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004}},{{0.703,0.002},{0.714,0.002},{0.774,0.004}}});
const vector<vector<ave_err_t>> Zt_ae({{{0.711,0.005},{0.724,0.004},{0.762,0.004}},{{0.700,0.003},{0.711,0.002},{0.767,0.002}}});
const int Za_seed[nbeta]={13124,862464,76753};
const int Zt_seed[nbeta]={5634,917453,324338};
dbvec_t Za,Zt;
boot_init_t bi;

dbvec_t alist(nbeta),zlist(nbeta);

//! hold the data for a single ensemble
class ens_data_t
{
public:
  size_t iult; //< input in the ultimate file
  size_t ib,T,L;
  int use_for_L;
  double aml;
  string path;
  
  size_t tmin[3],tmax[3];
  djack_t deltam_cr;
};
vector<ens_data_t> ens_data;
size_t nens_used;

//! read a single vector, for a specific mass and r, real or imaginary
inline djvec_t read(const char *what,const ens_data_t &ens,size_t im,size_t r,size_t reim)
{return read_djvec(combine("%s/data/corr%s",ens.path.c_str(),what),ens.T,ind({im,im,r,reim}));}

//! read a combination of r and return appropriately simmetrized
inline djvec_t read(const char *what,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  djvec_t o(ens.T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,ens,im,r,reim)*(r==0?1:rpar);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read averaging the three channels
inline djvec_t read(const char *what,const char *pat,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  return
    djvec_t(read(combine("%s_%c1%c1",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_%c2%c2",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	    read(combine("%s_%c3%c3",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim))
    /3.0;
}

//! read PP
inline djvec_t read_PP(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return read(combine("%s_%s",what,"P5P5").c_str(),ens,1,im,rpar,reim);}

//! read VV
inline djvec_t read_VV(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return read(what,"VV",ens,1,im,rpar,reim);}

//! read TV
inline djvec_t read_TV(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return -read(what,"TV",ens,-1,im,rpar,reim);}

//! read VT
inline djvec_t read_VT(const char *what,const ens_data_t &ens,size_t im,int rpar,size_t reim)
{return read(what,"VT",ens,-1,im,rpar,reim);}

//! compute the critical deltam
inline djack_t compute_deltam_cr(const ens_data_t &ens,size_t iq)
{
  djvec_t V0P5_LL=read("LL_V0P5",ens,-1,iq,-1,IM);
  djvec_t V0P5_0M=read("0M_V0P5",ens,-1,iq,-1,IM);
  djvec_t V0P5_0T=read("0T_V0P5",ens,-1,iq,-1,IM);
  djvec_t num_deltam_cr=forward_derivative(djvec_t(V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write(combine("%s/plots/num_deltam_cr.xmg",ens.path.c_str()));
  
  djvec_t V0P5_0P=read("0P_V0P5",ens,-1,iq,+1,RE);
  djvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write(combine("%s/plots/den_deltam_cr.xmg",ens.path.c_str()));
  
  djack_t deltam_cr=constant_fit(djvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),ens.tmin[iq],ens.tmax[iq],combine("%s/plots/deltam_cr_t.xmg",ens.path.c_str()));
  
  return deltam_cr;
}

//! read QED corrections
inline djvec_t read_QED(const char *pat,const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO)
{
  djvec_t c_0T=read("0T",pat,ens,tpar,im,1,RE);
  djvec_t c_0M=read("0M",pat,ens,tpar,im,1,RE);
  djvec_t c_LL=read("LL",pat,ens,tpar,im,1,RE);
  djvec_t(c_0T/c_LO).ave_err().write(combine("%s/plots/%s_0T.xmg",ens.path.c_str(),pat));
  djvec_t(c_0M/c_LO).ave_err().write(combine("%s/plots/%s_0M.xmg",ens.path.c_str(),pat));
  djvec_t(c_LL/c_LO).ave_err().write(combine("%s/plots/%s_LL.xmg",ens.path.c_str(),pat));
  djvec_t c=djvec_t(c_LL+2.0*djvec_t(c_0T+c_0M));
  
  djvec_t c_0P=read("0P",pat,ens,tpar,im,-1,IM);
  djvec_t(c_0P/c_LO).ave_err().write(combine("%s/plots/%s_0P.xmg",ens.path.c_str(),pat));
  djvec_t d=-(deltam_cr*c_0P);
  return djvec_t(c+2.0*d)*e2*sqr(eq[im]);
}

//! read for VV case
inline djvec_t read_QED_VV(const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO)
{return read_QED("VV",ens,tpar,im,deltam_cr,c_LO);}

//! read for TV case
inline djvec_t read_QED_TV(const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO)
{return read_QED("VV",ens,tpar,im,deltam_cr,c_LO);}

//! initialize gm2 calculation
void gm2_initialize(int narg,char **arg)
{
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  //read where to read input and how many ensemble
  string ens_pars=input.read<string>("UltimatePath");
  nm=input.read<size_t>("NMass");
  nr=input.read<size_t>("NR");
  ind.set_ranges({nm,nm,nr,2});
  init_common_IB(ens_pars);
  nens_used=input.read<int>("NEnsemble");
  Za.resize(nbeta);
  Zt.resize(nbeta);
  
  input.expect({"Ens","beta","L","UseForL","T","aml","tint_cr","tint_ss","tint_cc","path"});
  ens_data.resize(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      
      input.read(ens.iult);
      input.read(ens.ib);
      input.read(ens.L);
      input.read(ens.use_for_L);
      input.read(ens.T);
      input.read(ens.aml);
      for(size_t iq=0;iq<3;iq++)
	{
	  input.read(ens.tmin[iq]);
	  input.read(ens.tmax[iq]);
	}
      input.read(ens.path);
    }
}

//prepare the list of a and z
inline void prepare_az(int input_an_id)
{
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      const int imet=input_an_id/4;
      alist[ibeta]=1.0/lat_par[input_an_id].ainv[ibeta];
      zlist[ibeta]=lat_par[input_an_id].Z[ibeta];
      Za[ibeta].fill_gauss(Za_ae[imet][ibeta],Za_seed[ibeta]);
      Zt[ibeta].fill_gauss(Zt_ae[imet][ibeta],Zt_seed[ibeta]);
    }
}

#endif
