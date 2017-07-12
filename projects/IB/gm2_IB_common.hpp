#ifndef _GM2_IB_COMMON_HPP
#define _GM2_IB_COMMON_HPP

#include <common.hpp>
#include <gm2_IB_integrators.hpp>

const double M_P_phys[3]={0.775,0.686,2.9834};
const double M_V_phys[3]={0.775,1.0195,3.0969};
int use_extra_sources;
size_t nm,im,nr,include_ZA_perturb;
index_t ind_base,ind_extra;

//! systematic
const int nint_num_variations=4;
const int nfit_range_variations=2;
const int ncont_extrap=2;
const int nchir_variations=2;
const int nFSE_variations=3;

index_t ind_syst({
		  {"Input",ninput_an},
		  {"Chir",nchir_variations},
		  {"FSE",nFSE_variations},
		  {"Cont",ncont_extrap},
		  {"FitRn",nfit_range_variations},
		  {"IntNu",nint_num_variations}});
enum syst{c_input,c_chir,c_FSE,c_cont,c_fit_range,c_int_num};
template <syst comp> int case_of(int isyst){return ind_syst(isyst)[comp];}

//! hold the data for a single ensemble
class ens_data_t
{
public:
  size_t iult; //< input in the ultimate file
  size_t ib,T,L;
  double kappa;
  int use_for_L;
  double aml;
  string path;
  
  size_t tmin[3],tmax[3];
};
vector<ens_data_t> ens_data;
size_t nens_used;

//! return the bare mass given the ensemble and the quark index
inline double get_amq(const ens_data_t &ens,size_t im)
{
  const double ams[nbeta]={0.02363,0.02094,0.01612};
  const double amc[nbeta]={0.27903,0.24725,0.19037};
  double amq[3];
  amq[0]=ens.aml;
  amq[1]=ams[ens.ib];
  amq[2]=amc[ens.ib];
  
  return amq[im];
}

//! read the additional factor
size_t get_add_fact(const string &path,size_t im)
{
  static map<string,array<size_t,3>> facts;
  
  auto it=facts.emplace(path,array<size_t,3>{});
  if(it.second==true)
    {
      cout<<"Reading "<<endl;
      raw_file_t fin(path+"/data/add_fact.txt","r");
      for(size_t jm=0;jm<3;jm++) fin.read(it.first->second[jm]);
    }
  
  return it.first->second[im];
}

//! read a single vector, for a specific mass and r, real or imaginary
inline djvec_t read(const char *what,const ens_data_t &ens,size_t im,size_t r,size_t reim)
{
  string path=combine("%s/data/corr%s",ens.path.c_str(),what);
  string path_extra=path+"_"+qname[im][0]+qname[im][0];
  
  djvec_t out;
  
  if(file_exists(path))
    {
      out=read_djvec(path,ens.T,ind_base({im,im,r,reim}));
      if(use_extra_sources and file_exists(path_extra))
	{
	  size_t add_fact=get_add_fact(ens.path,im);
	  cout<<"Improving with "<<path_extra<<", weight: "<<add_fact<<endl;
	  
	  out+=add_fact*read_djvec(path_extra,ens.T,ind_extra({0,0,r,reim}));
	  out/=1+add_fact;
	}
    }
  else out=read_djvec(path_extra,ens.T,ind_extra({0,0,r,reim}));
  
  return out;
}

//! read a combination of r and return appropriately symmetrized
inline djvec_t read(const char *what,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  djvec_t o(ens.T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,ens,im,r,reim)*((r==0)?1:rpar);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read averaging the three channels
inline djvec_t read(const char *what,const char *pat,const ens_data_t &ens,int tpar,size_t im,int rpar,size_t reim)
{
  if(string(pat)=="VV")
    return
      djvec_t(read(combine("%s_%c1%c1",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	      read(combine("%s_%c2%c2",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim)+
	      read(combine("%s_%c3%c3",what,pat[0],pat[1]).c_str(),ens,tpar,im,rpar,reim))
      /3.0;
  else return read(combine("%s_%s",what,pat).c_str(),ens,tpar,im,rpar,reim);
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
  string ens_qpath=ens.path+"/plots_"+qname[iq];
  
  //compute M
  djvec_t P5P5_00=read("00_P5P5",ens,+1,iq,1,RE);
  djack_t M=constant_fit(effective_mass(P5P5_00),ens.tmin[iq],ens.tmax[iq],ens_qpath+"/P5P5_00_eff.xmg");
  cout<<qname[iq]<<" mass ens "<<ens.path<<": "<<M.ave_err()<<endl;
  
  //test jacknife per jacknife P5P5_00
  grace_file_t out(ens_qpath+"/P5P5_00_eff_jacks.xmg");
  for(size_t ijack=0;ijack<njacks;ijack++)
    {
      for(size_t t=0;t<P5P5_00.size();t++) out.write_xy(t,P5P5_00[t][njacks]*njacks-P5P5_00[t][ijack]*(njacks-1));
      out.new_data_set();
    }
  
  //measure mcrit according to eq.3 of hep-lat/0701012
  djvec_t V0P5_00=read("00_V0P5",ens,-1,iq,-1,IM);
  V0P5_00.ave_err().write(ens_qpath+"/V0P5_00.xmg");
  djvec_t m_cr_corr=forward_derivative(V0P5_00)/(2.0*P5P5_00);
  djack_t m_cr=constant_fit(m_cr_corr,ens.tmin[iq],ens.tmax[iq],ens_qpath+"/m_cr.xmg");
  effective_mass(V0P5_00,ens.T/2,-1).ave_err().write(ens_qpath+"/V0P5_00_eff.xmg");
  
  //load corrections
  djvec_t V0P5_LL=read("LL_V0P5",ens,-1,iq,-1,IM);
  djvec_t V0P5_0M=read("0M_V0P5",ens,-1,iq,-1,IM);
  djvec_t V0P5_0T=read("0T_V0P5",ens,-1,iq,-1,IM);
  djvec_t(V0P5_0M+V0P5_0T).ave_err().write(ens_qpath+"/V0P5_0MpT_norat.xmg");
  
  //build numerator
  djvec_t num_deltam_cr_corr=V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T);
  effective_mass(num_deltam_cr_corr,ens.T/2,-1).ave_err().write(ens_qpath+"/num_deltam_cr_corr_eff.xmg");
  djvec_t num_deltam_cr_alt_corr=num_deltam_cr_corr;
  num_deltam_cr_alt_corr.ave_err().write(ens_qpath+"/num_deltam_cr_alt_corr.xmg");
  djvec_t num_deltam_cr=forward_derivative(num_deltam_cr_corr);
  num_deltam_cr.ave_err().write(ens_qpath+"/num_deltam_cr.xmg");
  
  //load the derivative wrt counterterm
  djvec_t V0P5_0P=read("0P_V0P5",ens,-1,iq,+1,RE);
  //build denominator
  djvec_t den_deltam_cr_corr=V0P5_0P;
  effective_mass(den_deltam_cr_corr,ens.T/2,-1).ave_err().write(ens_qpath+"/den_deltam_cr_corr_eff.xmg");
  djvec_t den_deltam_cr=forward_derivative(den_deltam_cr_corr);
  den_deltam_cr.ave_err().write(ens_qpath+"/den_deltam_cr.xmg");
  
  //determine the counteterm
  djvec_t deltam_cr_corr=-num_deltam_cr/(2.0*den_deltam_cr);
  djack_t deltam_cr=constant_fit(deltam_cr_corr,ens.tmin[iq]-3,ens.tmax[iq],ens_qpath+"/deltam_cr_t.xmg");
  djvec_t deltam_cr_alt_corr=-num_deltam_cr_alt_corr/(2.0*den_deltam_cr_corr);
  djack_t deltam_cr_alt=constant_fit(deltam_cr_alt_corr,ens.tmin[iq]-6,ens.tmax[iq],ens_qpath+"/deltam_cr_alt_t.xmg");
  djack_t deltam_cr_alt_lin=poly_fit(deltam_cr_alt_corr,1,ens.tmin[iq]-3,ens.tmax[iq],ens_qpath+"/deltam_cr_alt_t_linfit.xmg")[0];
  
  // //check if there is a scalar missing
  // djvec_t V0P5_0S=read("0S_V0P5",ens,-1,iq,-1,IM);
  // forward_derivative(deltam_cr_alt_corr).ave_err().write(ens_qpath+"/deltam_cr_alt_der_t.xmg");
  // djvec_t deltam_cr_scal_corr=V0P5_0S/den_deltam_cr_corr;
  // deltam_cr_scal_corr.ave_err().write(ens_qpath+"/deltam_cr_scal_corr.xmg");
  // forward_derivative(deltam_cr_scal_corr).ave_err().write(ens_qpath+"/deltam_cr_scal_corr_der_t.xmg");
  // djvec_t deltam_cr_scal_coef_corr=forward_derivative(deltam_cr_alt_corr)/forward_derivative(deltam_cr_scal_corr);
  // djack_t deltam_cr_scal_coef=constant_fit(deltam_cr_scal_coef_corr,ens.tmin[iq],ens.tmax[iq],ens_qpath+"/deltam_cr_scal_contr_coef_t.xmg");
  // djvec_t deltam_cr_alt_corr_with_S=deltam_cr_alt_corr-deltam_cr_scal_coef*deltam_cr_scal_corr;
  // djack_t deltam_cr_alt_with_S=constant_fit(deltam_cr_alt_corr_with_S,ens.tmin[iq],ens.tmin[iq],ens_qpath+"/deltam_cr_alt_with_S_t.xmg");
  
  ///////////////////////// retuning of LO ///////////////////////////
  
  //how much should we retune kappa? (note that we should insert iP, so we put a minus, but this gets cancelled with the definition of dm_cr)
  djvec_t dm_cr_LO_corr=V0P5_00/(2.0*V0P5_0P);
  djack_t dm_cr_LO=constant_fit(dm_cr_LO_corr,ens.tmin[iq],ens.tmax[iq],ens_qpath+"/m_cr_retune.xmg");
  djack_t kappa_true=ens.kappa/(1+2*dm_cr_LO*ens.kappa);
  djack_t dk_cr_LO=kappa_true-ens.kappa;
  cout<<"True kappa cr (quark "<<qname[iq]<<") ens "<<ens_qpath<<": "<<smart_print(kappa_true.ave_err())
      <<" vs "<<ens.kappa<<", diff: "<<smart_print(dk_cr_LO.ave_err())<<endl;
  
  // //hep-lat/0101001 eq.2.5
  // double amq=get_amq(ens,iq);
  // djack_t tan_alpha_bare=amq/m_cr,arctan_alpha_bare=1/tan_alpha_bare;
  // djack_t alpha=atan2(amq,m_cr*Za_ae[0][ens.ib].ave())/M_PI*180;
  // cout<<"Twist angle in units of PI (quark "<<qname[iq]<<") ens "<<ens.path<<": "
  //     <<smart_print(alpha.ave_err())<<", am: "<<amq<<", mcr: "<<smart_print(m_cr.ave_err())<<endl;
  // cout<<"Contribution of scalar to pseudo (quark "<<qname[iq]<<") ens "<<ens.path<<": "<<smart_print(arctan_alpha_bare.ave_err())<<endl;
  // effective_mass(djvec_t(num_deltam_cr_corr+0*arctan_alpha_bare*V0P5_0S),ens.T/2,-1).ave_err().write(ens_qpath+"/num_deltam_cr_with_S_eff.xmg");
  
  //how much is the pion mass after retuning?
  djvec_t P5P5_0P=read("0P_P5P5",ens,+1,iq,-1,IM);
  djvec_t P5P5_00_retuned=P5P5_00+djack_t(2.0*dm_cr_LO)*P5P5_0P;
  djack_t M_retuned=constant_fit(effective_mass(P5P5_00_retuned),ens.tmin[iq],ens.tmax[iq],ens_qpath+"/P5P5_00_retuned_eff.xmg");
  cout<<"Variation of M due to k retuning: "<<smart_print(djack_t(M_retuned/M-1).ave_err())<<endl;
  
  return deltam_cr_alt;
}

//! read QED corrections
djack_t *Zm_fact;
inline djvec_t read_QED(const char *pat,const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO,const double a)
{
  string ens_qpath=ens.path+"/plots_"+qname[im];
  
  djvec_t c_0T=read("0T",pat,ens,tpar,im,1,RE);
  djvec_t c_0M=read("0M",pat,ens,tpar,im,1,RE);
  djvec_t c_LL=read("LL",pat,ens,tpar,im,1,RE);
  djvec_t(c_0T/c_LO).ave_err().write(ens_qpath+"/"+pat+"_0T.xmg");
  djvec_t(c_0M/c_LO).ave_err().write(ens_qpath+"/"+pat+"_0M.xmg");
  djvec_t(c_LL/c_LO).ave_err().write(ens_qpath+"/"+pat+"_LL.xmg");
  djvec_t c=djvec_t(c_LL+2.0*djvec_t(c_0T+c_0M));
  
  djvec_t c_0P=read("0P",pat,ens,tpar,im,-1,IM);
  c_0P.ave_err().write(ens_qpath+"/"+pat+"_0P_norat.xmg");
  djvec_t(c_0P/c_LO).ave_err().write(ens_qpath+"/"+pat+"_0P.xmg");
  djvec_t d=-(deltam_cr*c_0P);
  
  djack_t dmcrit_bare=deltam_cr*e2*sqr(eq[im]);
  djack_t kappa_new=ens.kappa/(1+2*dmcrit_bare*ens.kappa);
  cout<<"dmcrit_bare ("<<qname[im]<<" ens: "<<ens.path<<"): "<<dmcrit_bare.ave_err()<<", kappa: "<<kappa_new.ave_err()<<" old: "<<ens.kappa<<endl;
  
  djvec_t c_0S=read("0S",pat,ens,tpar,im,1,RE);
  c_0S.ave_err().write(ens_qpath+"/"+pat+"_0S_norat.xmg");
  djvec_t(c_0S/c_LO).ave_err().write(ens_qpath+"/"+pat+"_0S.xmg");
  double amq=get_amq(ens,im);
  
  double dm_bare_noe2=amq*(6.0*log(mu_MS*a)-22.596)/(16.0*sqr(M_PI));
  djvec_t e=(-dm_bare_noe2)*c_0S*(*Zm_fact); //minus? CHECK
  
  double dm_bare=dm_bare_noe2*e2*sqr(eq[im]);
  cout<<"dm_bare ("<<qname[im]<<" ens: "<<ens.path<<"): "<<dm_bare<<", total: "<<get_amq(ens,im)+dm_bare<<" old: "<<get_amq(ens,im)<<endl;
  
  return djvec_t(c+2.0*d+2.0*e)*e2*sqr(eq[im]);
}

//! read for PP case
inline djvec_t read_QED_PP(const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO,const double a)
{return read_QED("P5P5",ens,tpar,im,deltam_cr,c_LO,a);}

//! read for VV case
inline djvec_t read_QED_VV(const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO,const double a)
{return read_QED("VV",ens,tpar,im,deltam_cr,c_LO,a);}

//! read for TV case
inline djvec_t read_QED_TV(const ens_data_t &ens,const int tpar,const int im,const djack_t &deltam_cr,const djvec_t &c_LO,const double a)
{return read_QED("VV",ens,tpar,im,deltam_cr,c_LO,a);}

//! initialize gm2 calculation
inline void gm2_initialize(int narg,char **arg)
{
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  //read where to read input and how many ensemble
  string ens_pars=input.read<string>("UltimatePath");
  use_extra_sources=input.read<int>("UseExtraSources");
  nm=input.read<size_t>("NMass");
  im=input.read<size_t>("IMass");
  nr=input.read<size_t>("NR");
  include_ZA_perturb=input.read<size_t>("IncludeZAcorr");
  
  ind_base.set_ranges ({{"NMass",nm},{"NMass",nm},{"Nr",nr},{"RI",2}});
  ind_extra.set_ranges({{"NMass", 1},{"NMass", 1},{"Nr",nr},{"RI",2}});
  init_common_IB(ens_pars);
  nens_used=input.read<int>("NEnsemble");
  
  input.expect({"Ens","beta","L","UseForL","T","kappa","aml","tint_cr","tint_ss","tint_cc","path"});
  ens_data.resize(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      
      input.read(ens.iult);
      input.read(ens.ib);
      input.read(ens.L);
      input.read(ens.use_for_L);
      input.read(ens.T);
      input.read(ens.kappa);
      input.read(ens.aml);
      for(size_t iq=0;iq<3;iq++)
	{
	  input.read(ens.tmin[iq]);
	  input.read(ens.tmax[iq]);
	}
      input.read(ens.path);
    }
}
#endif
