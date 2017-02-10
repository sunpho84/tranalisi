#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <set>

#include <common.hpp>
#include <fit.hpp>
#include <functional>
#include <Kl2_IB_fit.hpp>

using namespace placeholders;
 
namespace
{
  const double a_cont=1.0e-5;
  const double inf_vol=1.0e10;
  
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

template <class TL,class Ta,class TD> TD FSE_dep_L4(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_savage(const TD &D1,const TD &D2, const TD &D3,const TD &M,double charge_radius,double pol_savage,const Ta &a,const TL &L)
{return e2*pow(charge_radius,2.0)/(L*a*L*a*L*a)/3.0*(D1*M+D2*4*M_PI*c_minus1_savage/(L*a))+D3*8*pow(M_PI,2.0)*M*c_minus1_savage*pol_savage/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class Tml,class TD> TD FSE_dep_pion_savage(const TD &D1,const TD &D2, const TD &D3,const TD &M,const TD &f0,const Tml &ml,const TD &ml_phys,const Ta &a,const TL &L)
{return e2*2/pow(4*M_PI*f0,2.0)*(l6-1-log(ml/ml_phys)+A_pion_charge_radius*(ml-ml_phys))/(L*a*L*a*L*a)/3.0*(D1*M+D2*4*M_PI*c_minus1_savage/(L*a))+D3*8*pow(M_PI,2.0)*M*c_minus1_savage*pion_pol_savage/(L*a*L*a*L*a*L*a);}


//////////////////////////////////////////////////// pion /////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &quad_dep,const Tml &ml,
			     const Tpars &ml_phys,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &L4dep,const Tpars &ML4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3+16*Cf04)*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_pion_savage(L3dep,L4dep,ML4dep,sqrtM2,f0,ml,ml_phys,a,L);
  else                 fitted_FSE=0;
  
  return e2*sqr(f0)*(4*Cf04+chir_dep+K*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{0.6,1.0},{0.0,10.0},boot_fit);
  pars.add_fsedep_pars({1.0,1.0},{1.0,1.0},{1.0,1.0},boot_fit);
  pars.add_LEC_pars({1.0e-5,1.0e-5},{6.0,1.0},{0.0,0.0},{0.0,10.0},{0.0,0.0},boot_fit);
  pars.add_ml_phys_par(ml_phys,boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  //boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,p[pars.iml_phys],ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iL4dep],p[pars.iML4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,pars.fit_ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2Pi<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		       pars.fit_ml_phys.ave(),pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.L4dep.ave(),pars.ML4dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,
		       pars.fit_ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				       without_with_fse*FSE_dep_pion_savage(pars.L3dep,pars.L4dep,pars.ML4dep,dboot_t(pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_f0,dboot_t(ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib]),pars.fit_ml_phys,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$(M^2_{\\pi^+}-M^2_{\\pi^0}) [GeV^2]");
  
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
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);
  
  return (4/3+3/(4*Cf04))*(chir_dep+Kpi*M2/den+Kk)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,{0.0,1.0},{0.0,1.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{1.5,0.18},{1.0,10.0},{1.0,10.0},{0.0,0.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  //boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
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
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// kaon_QED ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep,const Tm &ml,
				const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &L4dep,const Tpars &ML4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    sqrtM2K=pow(M2K,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-8/3*M2Pi/den*log(M2Pi/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2Pi/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_savage(L3dep,L4dep,ML4dep,sqrtM2K,kaon_charge_radius,kaon_pol_savage,a,L);
  else                 fitted_FSE=0;
  
  return e2*sqr(f0)*4*Cf04*(Kk+chir_dep+Kpi*M2Pi/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{2.24,1.0},{25.0,40.0},boot_fit);
  pars.add_fsedep_pars({0.23,0.05},{0.0,10.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({1.0e-4,1.0e-5},{20.0,4.0},{-50.0,10.0},{1.0,10.0},{0.0,0.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  //boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iL4dep],p[pars.iML4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag);
  cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;

  for(int univ_full_sub=0;univ_full_sub<2;univ_full_sub++)
    plot_chir_fit(combine(path.c_str(),(univ_full_sub==0)?"_univ":"_full"),ext_data,pars,
		  [&pars,&ms_phys,an_flag]
		  (double x,size_t ib)
		  {return cont_chir_ansatz_dM2K_QED<double,double,double>
		      (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),x,
		       ms_phys.ave(),pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.L4dep.ave(),pars.ML4dep.ave(),pars.adep_ml.ave(),an_flag);},
		  bind(cont_chir_ansatz_dM2K_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,_1,
		       ms_phys.ave(),a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag),
		  [&ext_data,&pars,univ_full_sub]
		  (size_t idata,bool without_with_fse,size_t ib)
		  {if(univ_full_sub==0)
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));
		      }
		    else
		      {
			return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wfse)/sqr(pars.fit_a[ib])-
				without_with_fse*FSE_dep_savage(pars.L3dep,pars.L4dep,pars.ML4dep,dboot_t(pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)),kaon_charge_radius,kaon_pol_savage,pars.fit_a[ib],ext_data[idata].L));
		      }
		  },
		  ml_phys,phys_res,"$$(M^2_{K^+}-M^2_{K^0})^{QED} [GeV^2]");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic linear fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_linear_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml;
  else                 chir_dep=0;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=sqr(a)*(adep+adep_ml*ml);
  else                 disc_eff=adep_ml*ml*sqr(a);
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag)
{
  //set_printlevel(3);
 
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0.0,1.0},{0.0,1.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},boot_fit);
  pars.add_adep_pars({0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_fsedep_pars({0.0,1.0},{0.0,0.0},{0.0,0.0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2Pi);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_linear_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_linear_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag,with_without_FSE]
		(double x,size_t ib)
		{return cont_chir_linear_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag,with_without_FSE);},
		bind(cont_chir_linear_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,yaxis_label);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic quadratic fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_quad_ansatz(const Tpars &C0,const Tpars &C1,const Tpars &C2,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml+C2*sqr(ml);
  else                 chir_dep=C1*ml;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=sqr(a)*(adep+adep_ml*ml);
  else                 disc_eff=adep_ml*ml*sqr(a);
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_quad_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0.0,1.0},{0.0,1.0},{0.0,0.0},{0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_adep_pars({0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_fsedep_pars({0.0,1.0},{0.0,0.0},{0.0,0.0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_quad_ansatz(p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_quad_ansatz(pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag,with_without_FSE]
		(double x,size_t ib)
		{return cont_chir_quad_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag,with_without_FSE);},
		bind(cont_chir_quad_ansatz<dboot_t,double,double>,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,yaxis_label);
  
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
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return Kpi0/(4*Cf04)*M2/den*(1+chir_dep)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{1.0,1.0},{-15.0,1.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{10.0,1.0},{0.0,10.0},{0.0,10.0},{0.0,0.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  //boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
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
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));},
		ml_phys,phys_res,"$$\\varepsilon_\\Pi0");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_K0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_K0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kk0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
				  const Ta &a,const Tpars &adep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=(3+16*Cf04)*M2/(4*Cf04*den)*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  return Kk0/(4*Cf04)*(1+chir_dep+K*M2/den)+disc_eff;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{-1.0,1.0},{1.0,1.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{1.0,10.0},{1.0,10.0},{1.0,10.0},{0.0,0.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  //boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon_K0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKK],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_K0(pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,an_flag);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_K0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KK.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon_K0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KK,pars.KPi,pars.K2Pi,_1,
		     a_cont,pars.adep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib]));},
		ml_phys,phys_res,"$$\\varepsilon_\\K0");
  
  return phys_res;
}

///////////////////////////////////////////////////////// kaon_QCD ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QCD(const Tpars &f0,const Tpars &B0,const Tpars &K,const Tpars &quad_dep,const Tm &ml,
				const Ta &a,const Tpars &adep,double L,const Tpars &fvedep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    M2=2*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep_ml*sqr(a)*ml;

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=fvedep*M2/den*exp(-sqrtM2*L)/pow(sqrtM2*L,1.5);
  else                 fitted_FSE=0;

  return B0*(1+chir_dep+K*M2/den)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QCD(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{1.0,1.0},{0.0,0.0},boot_fit);
  pars.add_fsedep_pars({0.0,1.0},{0.0,0.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({0.0,0.0},{1.0,10.0},{0.0,0.0},{1.0,10.0},{0.0,0.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iC);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,1.0,1.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QCD(p[pars.if0],p[pars.iB0],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QCD(pars.fit_f0,pars.fit_B0,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;

  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2K_QCD<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_dM2K_QCD<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.KPi,pars.K2Pi,_1,
		     a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pars.fit_z[ib]/sqr(pars.fit_a[ib])-
				without_with_fse*pars.L3dep*dboot_t(pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.25)/pow(4*M_PI*pars.fit_f0,2.0)*exp(-pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)*ext_data[idata].L)/pow(ext_data[idata].L,1.5)));},
		ml_phys,phys_res,"$$(M^2_{K^+}-M^2_{K^0})^{QCD}/(-2*Deltamud)[GeV]");
  
  return phys_res;
}
