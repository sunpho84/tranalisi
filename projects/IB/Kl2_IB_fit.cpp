#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <set>

#include <common.hpp>
#include <fit.hpp>
#include <functional>
#include <Kl2_IB_fit.hpp>

index_t ind_an;
using namespace placeholders;
 
namespace
{
  const double hslashc=0.197327;
  const double c_minus1_savage=-0.266596;
  const double pion_charge_radius=0.672/hslashc;
  const double kaon_charge_radius=0.56/hslashc;
  const double pion_pol_savage=7.1e-05/pow(hslashc,3.0);
  const double kaon_pol_savage=0.0;
  const double pol_savage_pi0=-1.1e-4/pow(hslashc,3.0);
  const double pol_savage_k0=-1.0e-4/pow(hslashc,3.0);
  const double l6=14.6;
  const double A_pion_charge_radius=-50.0;
  
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  
  const double mk0=0.497614;
  const double mkp=0.493677;
}

//! compute fitted piece of FSE
template <class TL,class Ta,class TD> TD FSE_dep_L3(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_ML3(const TD &D,const TD &M,const Ta &a,const TL &L)
{return D*M/(L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_L4(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_L6(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_savage(const TD &D1,const TD &D2, const TD &D3,const TD &M,double charge_radius,double pol_savage,const Ta &a,const TL &L)
{return e2*pow(charge_radius,2.0)/(L*a*L*a*L*a)/3.0*(D1*M+D2*4.0*M_PI*c_minus1_savage/(L*a))+D3*8.0*pow(M_PI,2.0)*M*c_minus1_savage*pol_savage/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class Tml,class TD> TD FSE_dep_pion_savage(const TD &D1,const TD &D2, const TD &D3,const TD &M,const TD &f0,const Tml &ml,const TD &ml_phys,const Ta &a,const TL &L)
{return e2*2.0/pow(4.0*M_PI*f0,2.0)*(l6-1.0-log(ml/ml_phys)+A_pion_charge_radius*(ml-ml_phys))/(L*a*L*a*L*a)/3.0*(D1*M+D2*4.0*M_PI*c_minus1_savage/(L*a))+D3*8.0*pow(M_PI,2.0)*M*c_minus1_savage*pion_pol_savage/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class Tml,class TD> TD FSE_dep_pion_savage_new(const TD &D1,const TD &M,const TD &f0,const Tml &ml,const TD &ml_phys,const Ta &a,const TL &L)
{return e2*2.0/pow(4.0*M_PI*f0,2.0)*(l6-1.0-log(ml/ml_phys)+A_pion_charge_radius*(ml-ml_phys))/(L*a*L*a*L*a)/3.0*D1*M;}

template <class TL,class Ta,class TD> TD FSE_dep_pion_savage_new_exp(const TD &D1,const TD &M,const Ta &a,const TL &L)
{return e2*pow(pion_charge_radius,2.0)/(L*a*L*a*L*a)/3.0*D1*M;}

template <class TL,class Ta,class TD> TD FSE_dep_kaon_savage(const TD &D1,const TD &M,const Ta &a,const TL &L)
{return e2*pow(kaon_charge_radius,2.0)/(L*a*L*a*L*a)/3.0*D1*M;}

//////////////////////////////////////////////////// pion /////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &quad_dep,const Tml &ml,
			     const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3.0+16.0*Cf04)*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage_new_exp(L3dep,sqrtM2,a,L);
  else                 fitted_FSE=0.0;
  
  return e2*sqr(f0)*(4.0*Cf04+chir_dep+K*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  
  ave_err_t adep_Guess={1.2e-3,1.0e-3};
  ave_err_t adep_ml_Guess={0.1,0.01};
  ave_err_t C_Guess={4.0e-5,2.5e-6};
  ave_err_t KPi_Guess={6.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={1.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2Pi<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),FSE_an(an_flag));},
		  bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_pion_savage_new_exp(pars.L3dep,dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$M^2_{\\pi^+}-M^2_{\\pi^0} [GeV^2]",beta_list,univ_full_sub,FSE_an(an_flag));
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// kaon_QED ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
				const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2.0*B0*ml,
    M2K=B0*(ml+ms),
    sqrtM2K=pow(M2K,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-8.0/3.0*M2Pi/den*log(M2Pi/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2Pi/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_kaon_savage(L3dep,sqrtM2K,a,L);
  else                 fitted_FSE=0.0;
  
  return e2*sqr(f0)*4.0*Cf04*(Kk+chir_dep+Kpi*M2Pi/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={4.0e-3,2.0e-3};
  ave_err_t adep_ml_Guess={0.06,0.03};
  ave_err_t C_Guess={1.0e-4,4.0e-5};
  ave_err_t KPi_Guess={8.0,2.0};
  ave_err_t KK_Guess={0.6,0.2};
  ave_err_t K2Pi_Guess={-15.0,6.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={1.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,&ms_phys,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2K_QED<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		       ms_phys.ave(),pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2K_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       ms_phys.ave(),a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_kaon_savage(pars.L3dep,dboot_t(pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$[M^2_{K^+}-M^2_{K^0}]_{QED} [GeV^2]",beta_list,univ_full_sub,FSE_an(an_flag));
  
  return phys_res;
}

///////////////////////////////////////////////////////// kaon_QCD ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QCD(const Tpars &f0,const Tpars &B0,const Tpars &Kc,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
				const Ta &a,const Tpars &adep,double L,const Tpars &fvedep,const size_t an_flag)
{
  Tpars
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=fvedep*M2/den*exp(-sqrtM2*L)/pow(sqrtM2*L,1.5);
  else                 fitted_FSE=0.0;

  return Kc*(1.0+chir_dep)+K*M2/den+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QCD(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={1.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={0.0,0.1};
  ave_err_t KPi_Guess={1.0,1.0};
  ave_err_t KK_Guess={2.0,0.1};
  ave_err_t K2Pi_Guess={1.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iC);
  //boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,1.0,1.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QCD(p[pars.if0],p[pars.iB0],p[pars.iKK],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QCD(pars.fit_f0,pars.fit_B0,pars.KK,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag);
  cout<<"dM2K_QCD_over_minus_two_Deltamud: "<<phys_res.ave_err()<<endl;

  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2K_QCD<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.KK.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_dM2K_QCD<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.KK,pars.KPi,pars.K2Pi,_1,
		     a_cont,pars.adep,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pars.fit_z[ib]/pars.fit_a[ib]-
				without_with_fse*pars.L3dep*dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.25)/pow(4.0*M_PI*pars.fit_f0,2.0)*exp(-pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)*ext_data[idata].L)/pow(ext_data[idata].L,1.5)));},
		ml_phys,phys_res,"$$[M^2_{K^+}-M^2_{K^0}]_{QCD}/(-2\\Delta m_{ud}) [GeV]",beta_list);
  
  return phys_res;
}

///////////////////////////////////////////////////////// M2Pi0g ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_M2Pi0g(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
			      const Ta &a,const Tpars &adep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*M2/den;

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return Kpi0/(4.0*Cf04)*M2/den*(1.0+chir_dep)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_M2Pi0g(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={1.5e-4,1.0e-4};
  ave_err_t adep_ml_Guess={0.015,0.01};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={6.2e-3,4.0e-3};
  ave_err_t KK_Guess={0.0,1.0};
  ave_err_t K2Pi_Guess={0.8,0.5};
  ave_err_t K2K_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKK,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_M2Pi0g(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_M2Pi0g(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,an_flag);
  cout<<"M2Pi0g: "<<phys_res.ave_err()<<endl;

  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_M2Pi0g<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_M2Pi0g<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		     a_cont,pars.adep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));},
		ml_phys,phys_res,"$$[\\delta M^2_{\\pi^0}]_{QED} [GeV^2]",beta_list);
  
  return phys_res;
}

///////////////////////////////////////////////////////// M2K0g ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_M2K0g(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kk0,const Tpars &K,const Tpars &Klog,const Tm &ml,
			     const Ta &a,const Tpars &adep,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=Klog*(3.0+16.0*Cf04)*M2/(4.0*Cf04*den)*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=0.0;

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;

  return Kk0/(4.0*Cf04)*(1.0+chir_dep+K*M2/den)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_M2K0g(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={4.0e-4,2.0e-4};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={9.0e-5,1.0e-5};
  ave_err_t KPi_Guess={2.6,1.0};
  ave_err_t KK_Guess={4.0e-4,2.0e-4};
  ave_err_t K2Pi_Guess={0.0,0.1};
  ave_err_t K2K_Guess={-0.3,0.2};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iK2K,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2Pi);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_M2K0g(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKK],p[pars.iKPi],p[pars.iK2K],ml,ac,p[pars.iadep],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_M2K0g(pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2K,ml_phys,a_cont,pars.adep,an_flag);
  cout<<"M2K0g: "<<phys_res.ave_err()<<endl;

  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_M2K0g<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KK.ave(),pars.KPi.ave(),pars.K2K.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),an_flag);},
		bind(cont_chir_ansatz_M2K0g<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2K,_1,
		     a_cont,pars.adep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));},
		ml_phys,phys_res,"$$[\\delta M^2_{K^0}]_{QED} [GeV^2]",beta_list);
  
  return phys_res;
}
/*
///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
			       const Ta &a,const Tpars &adep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return (4.0/3.0+3.0/(4.0*Cf04))*(Kk+chir_dep+Kpi*M2/den)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={-0.3,0.15};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={1.4,0.18};
  ave_err_t KK_Guess={0.2,0.05};
  ave_err_t K2Pi_Guess={8.0,1.5};
  ave_err_t K2K_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,an_flag);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;

  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		     a_cont,pars.adep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma",beta_list);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_Pi0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_Pi0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
				   const Ta &a,const Tpars &adep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*M2/den;

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return Kpi0/(4.0*Cf04)*M2/den*(1.0+chir_dep)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.3,0.1};
  ave_err_t adep_ml_Guess={-35.0,8.0};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={9.0,2.0};
  ave_err_t KK_Guess={0.8,0.2};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKK,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon_Pi0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_Pi0(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,an_flag);
  cout<<"epsilon_Pi0: "<<phys_res.ave_err()<<", exp: 0.07"<<endl;

 
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_Pi0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon_Pi0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		     a_cont,pars.adep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\pi^0",beta_list);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_K0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_K0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kk0,const Tpars &K,const Tpars &quad_dep,const Tpars &Klog,const Tm &ml,
				  const Ta &a,const Tpars &adep,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=Klog*(3.0+16.0*Cf04)*M2/(4.0*Cf04*den)*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;

  return Kk0/(4.0*Cf04)*(1.0+chir_dep+K*M2/den)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={17.6,0.2};
  ave_err_t KK_Guess={0.5,0.1};
  ave_err_t K2Pi_Guess={21.0,6.0};
  ave_err_t K2K_Guess={1.0,1.0};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iK2K,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iC);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon_K0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKK],p[pars.iKPi],p[pars.iK2Pi],p[pars.iK2K],ml,ac,p[pars.iadep],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_K0(pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,pars.K2K,ml_phys,a_cont,pars.adep,an_flag);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;

  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_K0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KK.ave(),pars.KPi.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon_K0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,pars.K2K,_1,
		     a_cont,pars.adep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\K^0",beta_list);
  
  return phys_res;
}
*/
///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
			       const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage_new_exp(L3dep,sqrtM2,a,L);
  else                 fitted_FSE=0.0;

  return (4.0/3.0+3.0/(4.0*Cf04))*(Kk+chir_dep+Kpi*M2/den)*(1.0-fitted_FSE)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={-0.3,0.15};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={1.4,0.18};
  ave_err_t KK_Guess={0.2,0.05};
  ave_err_t K2Pi_Guess={8.0,1.5};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={1.0,10.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*
				       (1.0+without_with_fse*FSE_dep_pion_savage_new_exp(pars.L3dep,dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L)));
		      }
		  },
		  ml_phys,phys_res,"$$\\varepsilon_\\gamma",beta_list,univ_full_sub,an_flag);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_Pi0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_Pi0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
				   const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*M2/den;

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage_new_exp(L3dep,sqrtM2,a,L);
  else                 fitted_FSE=0.0;
  
  return Kpi0/(4.0*Cf04)*M2/den*(1.0+chir_dep)*(1.0-fitted_FSE)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.3,0.1};
  ave_err_t adep_ml_Guess={-35.0,8.0};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={9.0,2.0};
  ave_err_t KK_Guess={0.8,0.2};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={-40.0,8.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKK,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon_Pi0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_Pi0(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"epsilon_Pi0: "<<phys_res.ave_err()<<", exp: 0.07"<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_epsilon_Pi0<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_epsilon_Pi0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*
				       (1.0+without_with_fse*FSE_dep_pion_savage_new_exp(pars.L3dep,dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L)));
		      }
		  },
		  ml_phys,phys_res,"$$\\varepsilon_\\pi^0",beta_list,univ_full_sub,an_flag);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_K0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_K0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kk0,const Tpars &K,const Tpars &quad_dep,const Tpars &Klog,const Tm &ml,
				  const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=Klog*(3.0+16.0*Cf04)*M2/(4.0*Cf04*den)*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage_new_exp(L3dep,sqrtM2,a,L);
  else                 fitted_FSE=0.0;
  
  return Kk0/(4.0*Cf04)*(1.0+chir_dep+K*M2/den)*(1.0-fitted_FSE)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={17.6,0.2};
  ave_err_t KK_Guess={0.5,0.1};
  ave_err_t K2Pi_Guess={21.0,6.0};
  ave_err_t K2K_Guess={1.0,1.0};
  ave_err_t L3dep_Guess={-60.0,10.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iK2K,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon_K0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKK],p[pars.iKPi],p[pars.iK2Pi],p[pars.iK2K],ml,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_K0(pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,pars.K2K,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_epsilon_K0<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KK.ave(),pars.KPi.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		  bind(cont_chir_ansatz_epsilon_K0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,pars.K2K,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*
				       (1.0+without_with_fse*FSE_dep_pion_savage_new_exp(pars.L3dep,dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L)));
		      }
		  },
		  ml_phys,phys_res,"$$\\varepsilon_\\K^0",beta_list,univ_full_sub,an_flag);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic linear fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class TM,class Ta>
Tpars cont_chir_linear_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const TM &MD,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml;
  else                 chir_dep=0.0;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_ML3(L3dep,MD,a,L);
      else                 fitted_FSE=0.0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
 
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={0.0,1.0};
  ave_err_t KPi_Guess={0.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,0.1};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_az_pars(a,z,boot_fit);
  pars.add_adep_pars(adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKPi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(with_without_FSE==0)
    {
      boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  else
    {
      if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2Pi);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_linear_ansatz(p[pars.iC],p[pars.iKPi],ml,MD,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag,with_without_FSE);},cov_flag);

  dboot_t phys_res=cont_chir_linear_ansatz(pars.C,pars.KPi,ml_phys,dboot_t(ext_data[0].aMaux/pars.fit_a[0]),a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE);

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,&ext_data,an_flag,with_without_FSE]
		  (double x,size_t ib)
		  {return cont_chir_linear_ansatz<double,double,double>
		      (pars.C.ave(),pars.KPi.ave(),x,dboot_t(ext_data[0].aMaux/pars.fit_a[0])[0],pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag,with_without_FSE);},		  
		  bind(cont_chir_linear_ansatz<dboot_t,double,dboot_t,double>,pars.C,pars.KPi,_1,dboot_t(ext_data[0].aMaux/pars.fit_a[0]),a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE),
		  [&ext_data,apow,zpow,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				       without_with_fse*FSE_dep_ML3(pars.L3dep,dboot_t(ext_data[idata].aMaux/pars.fit_a[ib]),pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,yaxis_label,beta_list,univ_full_sub);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic quadratic fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_quad_ansatz(const Tpars &C0,const Tpars &C1,const Tpars &C2,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml+C2*sqr(ml);
  else                 chir_dep=C1*ml;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0.0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_quad_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={0.0,1.0};
  ave_err_t KPi_Guess={0.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_az_pars(a,z,boot_fit);
  pars.add_adep_pars(adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(with_without_FSE==0)
    {
      boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  else
    {
      if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_quad_ansatz(p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_quad_ansatz(pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE);

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag,with_without_FSE]
		  (double x,size_t ib)
		  {return cont_chir_quad_ansatz<double,double,double>
		      (pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag,with_without_FSE);},
		  bind(cont_chir_quad_ansatz<dboot_t,double,double>,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE),
		  [&ext_data,apow,zpow,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				       without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,yaxis_label,beta_list,univ_full_sub);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic constant fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_constant_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml;
  else                 chir_dep=0.0;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0.0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_constant_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
 
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={0.0,1.0};
  ave_err_t KPi_Guess={0.0,0.1};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,0.1};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_az_pars(a,z,boot_fit);
  pars.add_adep_pars(adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  //if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKPi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(with_without_FSE==0)
    {
      boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  else
    {
      if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKPi);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2Pi);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_constant_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_constant_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE);

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag,with_without_FSE]
		  (double x,size_t ib)
		  {return cont_chir_constant_ansatz<double,double,double>
		      (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag,with_without_FSE);},
		  bind(cont_chir_constant_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE),
		  [&ext_data,apow,zpow,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				       without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,yaxis_label,beta_list,univ_full_sub);
  
  return phys_res;
}
/*
//////////////////////////////////////////////////// pion savage new exp/////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &quad_dep,const Tml &ml,
			     const Tpars &ml_phys,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3.0+16.0*Cf04)*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage_new(L3dep,sqrtM2,f0,ml,ml_phys,a,L);
  else                 fitted_FSE=0.0;
  
  return e2*sqr(f0)*(4.0*Cf04+chir_dep+K*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  
  ave_err_t adep_Guess={1.2e-3,1.0e-3};
  ave_err_t adep_ml_Guess={0.1,0.01};
  ave_err_t C_Guess={4.0e-5,2.5e-6};
  ave_err_t KPi_Guess={6.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={3.0,2.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_ml_phys_par(ml_phys,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,p[pars.iml_phys],ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,pars.fit_ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2Pi<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		       pars.fit_ml_phys.ave(),pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,
		       pars.fit_ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_pion_savage_new(pars.L3dep,dboot_t(pow(2.0*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_f0,dboot_t(ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib]),pars.fit_ml_phys,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$M^2_{\\pi^+}-M^2_{\\pi^0} [GeV^2]",beta_list,univ_full_sub,an_flag);
  
  return phys_res;
}

//////////////////////////////////////////////////// pion L6/////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &quad_dep,const Tml &ml,
			     const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3.0+16.0*Cf04)*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L6(L3dep,a,L);
  else                 fitted_FSE=0.0;
  
  return e2*sqr(f0)*(4.0*Cf04+chir_dep+K*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  
  ave_err_t adep_Guess={1.2e-3,1.0e-3};
  ave_err_t adep_ml_Guess={0.1,0.01};
  ave_err_t C_Guess={4.0e-5,2.5e-6};
  ave_err_t KPi_Guess={6.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,1.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={3.0,2.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2Pi<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_L6(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$M^2_{\\pi^+}-M^2_{\\pi^0} [GeV^2]",beta_list,univ_full_sub);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// kaon_QED L6////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
				const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-8.0/3.0*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L6(L3dep,a,L);
  else                 fitted_FSE=0.0;
  
  return e2*sqr(f0)*4.0*Cf04*(Kk+chir_dep+Kpi*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={4.0e-3,2.0e-3};
  ave_err_t adep_ml_Guess={0.06,0.03};
  ave_err_t C_Guess={1.0e-4,4.0e-5};
  ave_err_t KPi_Guess={8.0,2.0};
  ave_err_t KK_Guess={0.6,0.2};
  ave_err_t K2Pi_Guess={-15.0,6.0};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={1.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2K_QED<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2K_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_L6(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$[M^2_{K^+}-M^2_{K^0}]_{QED} [GeV^2]",beta_list,univ_full_sub);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
			       const Ta &a,const Tpars &adep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return (4.0/3.0+3.0/(4.0*Cf04))*(Kk+chir_dep+Kpi*M2/den)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={-0.3,0.15};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={1.4,0.18};
  ave_err_t KK_Guess={0.2,0.05};
  ave_err_t K2Pi_Guess={1.0,0.1};
  ave_err_t K2K_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,an_flag);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;

  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		     a_cont,pars.adep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma",beta_list);
  
  return phys_res;
}

///////////////////////////////////////////////////////// M2Pi0g ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_M2Pi0g(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
			      const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2.0*B0*ml,
    den=sqr(Tpars(4.0*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*M2/den;

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
  else                 fitted_FSE=0.0;
  
  return Kpi0/(4.0*Cf04)*M2/den*(1.0+chir_dep)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_M2Pi0g(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={1.5e-4,1.0e-4};
  ave_err_t adep_ml_Guess={0.015,0.01};
  ave_err_t C_Guess={5.0e-5,1.0e-5};
  ave_err_t KPi_Guess={6.2e-3,4.0e-3};
  ave_err_t KK_Guess={0.0,1.0};
  ave_err_t K2Pi_Guess={0.8,0.5};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==1) boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKK,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_ansatz_M2Pi0g(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_M2Pi0g(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"M2Pi0g: "<<phys_res.ave_err()<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_M2Pi0g<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		       pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_M2Pi0g<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$[\\delta M^2_{\\pi^0}]_{QED} [GeV^2]",beta_list,univ_full_sub);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic linear fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_linear_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml;
  else                 chir_dep=0.0;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a);
  else                 disc_eff=0.0;
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0.0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
 
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);

  ave_err_t adep_Guess={0.0,1.0};
  ave_err_t adep_ml_Guess={0.0,0.1};
  ave_err_t C_Guess={0.0,1.0};
  ave_err_t KPi_Guess={0.0,1.0};
  ave_err_t KK_Guess={0.0,0.1};
  ave_err_t K2Pi_Guess={0.0,0.1};
  ave_err_t K2K_Guess={0.0,0.1};
  ave_err_t L3dep_Guess={0.0,1.0};
  ave_err_t L4dep_Guess={0.0,0.1};
  ave_err_t ML4dep_Guess={0.0,0.1};
  
  pars.add_az_pars(a,z,boot_fit);
  pars.add_adep_pars(adep_Guess,adep_ml_Guess,boot_fit);
  pars.add_LEC_pars(C_Guess,KPi_Guess,KK_Guess,K2Pi_Guess,K2K_Guess,boot_fit);
  pars.add_fsedep_pars(L3dep_Guess,L4dep_Guess,ML4dep_Guess,boot_fit);
  if(chir_an(an_flag)==0) boot_fit.fix_par_to(pars.iKPi,0.0);
  if(cont_an(an_flag)==0) boot_fit.fix_par_to(pars.iadep,0.0);
  if(with_without_FSE==0)
    {
      boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  else
    {
      if(FSE_an(an_flag)==0)  boot_fit.fix_par_to(pars.iL3dep,0.0);
    }
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2Pi);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MD,double ac,double L)
			 {return cont_chir_linear_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],an_flag,with_without_FSE);},cov_flag);

  dboot_t phys_res=cont_chir_linear_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE);

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,&ext_data,an_flag,with_without_FSE]
		  (double x,size_t ib)
		  {return cont_chir_linear_ansatz<double,double,double>
		      (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),an_flag,with_without_FSE);},		  
		  bind(cont_chir_linear_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,an_flag,with_without_FSE),
		  [&ext_data,apow,zpow,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				       without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,yaxis_label,beta_list,univ_full_sub);
  
  return phys_res;
}
*/
