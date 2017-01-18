#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fit.hpp>
#include <functional>
#include <Kl2_IB_fit.hpp>

using namespace placeholders;

namespace
{
  const double a_cont=1e-5;
  const double inf_vol=1e10;
  
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  
  const double mk0=0.497614;
  const double mkp=0.493677;
}

//! compute fitted piece of FSE
template <class TL,class Ta,class TD> TD FSE_dep(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_sav(const TD &D,const TD &M,const Ta &a,const TL &L)
{return D*(M/(L*a*L*a*L*a)-4*M_PI*0.266596/(L*a*L*a*L*a*L*a)-0.266596*8*pow(M_PI,2.0)*M*7e-05/pow(0.19731,3.0)/(L*a*L*a*L*a*L*a)/(e2*pow(0.672/0.19731,2.0)/3));}

//! holds index and out pars
class cont_chir_fit_pars_t
{
public:
  dbvec_t ori_a,ori_z;
  dboot_t ori_f0,ori_B0;
  size_t nbeta;
  vector<size_t> ipara,iparz;
  size_t if0,iB0;
  size_t iadep,iadep_ml,iL3dep;
  size_t iC,iKPi,iKK;
  bool fitting_a,fitting_z;
  dbvec_t fit_a,fit_z;
  dboot_t fit_f0,fit_B0;
  dboot_t adep,adep_ml,L3dep;
  dboot_t C,KPi,KK;
  
  cont_chir_fit_pars_t(size_t nbeta)
    :
    nbeta(nbeta),
    ipara(nbeta),
    iparz(nbeta),
    fitting_a(true),
    fitting_z(true),
    fit_a(nbeta),
    fit_z(nbeta)
  {}
  
  //! return a
  double get_a(const vector<double> &p,size_t ib,size_t iel) const
  {
    if(fitting_a) return p[ipara[ib]];
    else          return ori_a[ib][iel];
  }
  
  //! return z
  double get_z(const vector<double> &p,size_t ib,size_t iel) const
  {
    if(fitting_z) return p[iparz[ib]];
    else          return ori_z[ib][iel];
  }
  
  //! add a and az to be self-fitted
  void add_az_pars(const dbvec_t &a,const dbvec_t &z,boot_fit_t &boot_fit)
  {
    ori_a=a;
    ori_z=z;
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      {
	ipara[ibeta]=boot_fit.add_self_fitted_point(fit_a[ibeta],combine("a[%zu]",ibeta),a[ibeta],DO_NOT_CORRELATE);
	iparz[ibeta]=boot_fit.add_self_fitted_point(fit_z[ibeta],combine("z[%zu]",ibeta),z[ibeta],DO_NOT_CORRELATE);
      }
  }
  
  //! add adep and adep_ml
  void add_adep_pars(const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,boot_fit_t &boot_fit)
  {
    // adep, adep_ml
    iadep=boot_fit.add_fit_par(adep,"adep",adep_guess.ave,adep_guess.err);
    iadep_ml=boot_fit.add_fit_par(adep_ml,"adep_ml",adep_ml_guess.ave,adep_ml_guess.err);
  }
  
  //! add L3dep
  void add_L3dep_par(const ave_err_t &L3dep_guess,boot_fit_t &boot_fit)
  {
    iL3dep=boot_fit.add_fit_par(L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  }
  
  //! add all common pars
  void add_common_pars(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
		       const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &L3dep_guess,boot_fit_t &boot_fit)
  {
    add_az_pars(a,z,boot_fit);
    add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
    //f0 and B0
    ori_f0=f0;
    ori_B0=B0;
    if0=boot_fit.add_self_fitted_point(fit_f0,"f0",f0,DO_NOT_CORRELATE);
    iB0=boot_fit.add_self_fitted_point(fit_B0,"B0",B0,DO_NOT_CORRELATE);
    
    add_L3dep_par(L3dep_guess,boot_fit);
  }
  
  //! add the low energy constants
  void add_LEC_pars(const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &KK_guess,boot_fit_t &boot_fit)
  {
    iC=boot_fit.add_fit_par(C,"C",C_guess.ave,C_guess.err);
    iKPi=boot_fit.add_fit_par(KPi,"KPi",KPi_guess.ave,KPi_guess.err);
    iKK=boot_fit.add_fit_par(KK,"KK",KK_guess.ave,KK_guess.err);
  }
  
  //////////////////////////////////////////////////////////////////////
  
  void print_common_pars() const
  {
    for(size_t ia=0;ia<ori_a.size();ia++)
      cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<ori_a[ia].ave_err()<<", ratio: "<<dboot_t(ori_a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
    for(size_t iz=0;iz<ori_z.size();iz++)
      cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<ori_z[iz].ave_err()<<", ratio: "<<dboot_t(ori_z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
    //
    cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<ori_f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/ori_f0-1).ave_err()<<endl;
    cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<ori_B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/ori_B0-1).ave_err()<<endl;
    //
    cout<<"Adep: "<<adep.ave_err()<<endl;
    cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
    cout<<"L3dep: "<<L3dep.ave_err()<<endl;
  }
  
  //! print the value
  void print_LEC_pars() const
  {
    cout<<"C: "<<C.ave_err()<<endl;
    cout<<"KPi: "<<KPi.ave_err()<<endl;
    cout<<"KK: "<<KK.ave_err()<<endl;
  }
};

//! plot the continuum-chiral extrapolation
void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label)
{
  //search max renormalized mass
  double ml_max=0;
  for(auto &data : ext_data)
    ml_max=max(ml_max,dboot_t(data.aml/pars.fit_a[data.ib]/pars.fit_z[data.ib]).ave());
  ml_max*=1.1;
  
  //prepare plot
  grace_file_t fit_file(path);
  fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<pars.fit_a.size();ib++)
	{
	  fit_file.new_data_set();
	  //put data without fse to brown
	  if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	  
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      fit_file<<dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave()<<" "<<
		fun_data(idata,without_with_fse,ib).ave_err()<<endl;
	}
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
}

//! perform the fit to the continuum limit
void cont_chir_fit_minimize
(const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,boot_fit_t &boot_fit,double apow,double zpow,
 const function<double(const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)> &cont_chir_ansatz,bool cov_flag)
{
  //set_printlevel(3);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    boot_fit.add_point(//numerical data
		       [&ext_data,&pars,idata,apow,zpow]
		       (vector<double> p,int iel) //dimension 2
		       {return ext_data[idata].wfse[iel]*pow(pars.get_z(p,ext_data[idata].ib,iel),zpow)/pow(pars.get_a(p,ext_data[idata].ib,iel),apow);},
		       //ansatz
		       [idata,&pars,&ext_data,&cont_chir_ansatz]
		       (vector<double> p,int iel)
		       {
			 size_t ib=ext_data[idata].ib;
			 double ac=pars.get_a(p,ib,iel);
			 double zc=pars.get_z(p,ib,iel);
			 double ml=ext_data[idata].aml/ac/zc;
			 double ms=ext_data[idata].ams/ac/zc;
			 double L=ext_data[idata].L;
			 return cont_chir_ansatz(p,pars,ml,ms,ac,L);
		       },
		       //for covariance/error
		       dboot_t(ext_data[idata].wfse/pow(pars.ori_a[ext_data[idata].ib],apow)),1/*correlate*/);
  
  //! fit
  boot_fit.fit(cov_flag);
  
  //print parameters
  pars.print_common_pars();
  pars.print_LEC_pars();
}

//////////////////////////////////////////////////// pion /////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tml &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    sqrtM2=pow(2*B0*ml,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=-(3+16*Cf04)*M2/den*log(M2);
  else        chir_log=0;
  
  Tpars fitted_FSE=e2*pow(0.672/0.19731,2.0)/3*FSE_dep_sav(L3dep,sqrtM2,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return e2*sqr(f0)*(4*Cf04+chir_log+K*M2/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{0.6,1},{0.0,10},{0.0,1},boot_fit);
  pars.add_LEC_pars({1e-5,1e-5},{6,1},{0,0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  //boot_fit.fix_par(pars.iadep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[chir_an](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
	 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,B0,pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2Pi<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib])-
				without_with_fse*e2*pow(0.672/0.19731,2.0)/3*FSE_dep_sav(pars.L3dep,dboot_t(pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,"$$(M^2_{\\pi^+}-M^2_{\\pi^0}) [GeV^2]");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,
			       const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  // if(chir_an)
  chir_log=M2Pi/den*log(M2Pi)-M2K/den*log(M2K);
  // else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return (2+3/(4*Cf04))*(chir_log+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,{0,1},{0.0,0.0},{-78,14},boot_fit);
  pars.add_LEC_pars({5e-5,1e-5},{1.5,0.18},{-1.2,0.16},boot_fit);
  boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[chir_an](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
     {return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);}
			 ,cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(f0,B0,pars.C,pars.KPi,pars.KK,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// kaon_QED ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,
			       const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=-(3+8*Cf04)*M2K/den*log(M2K)-8*Cf04*M2Pi/den*log(M2Pi);
  else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L)*pow(M2K,0.5);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return e2*sqr(f0)*(4*Cf04+chir_log+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{2.24,1},{25,40},{0.23,0.05},boot_fit);
  pars.add_LEC_pars({1e-4,1e-5},{20,4},{-50,10},boot_fit);
  boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[chir_an](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QED(f0,B0,pars.C,pars.KPi,pars.KK,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2K_QED<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_dM2K_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib])-
				without_with_fse*FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L)*pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5));},
		ml_phys,phys_res,"$$(M^2_{K^+}-M^2_{K^0})^{QED} [GeV^2]");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic linear fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_linear_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const Ta &a,const Tpars &adep,const Tpars &adep_ml)
{return C0+C1*ml+a*a*(adep+adep_ml*ml);}

//! perform the fit to the continuum limit
dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,bool fix_adep_ml,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0,1},{0,1},{0,0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  pars.add_adep_pars({0,1},{0,0},boot_fit);
  if(fix_adep_ml) boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_linear_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],p[pars.iadep_ml]);},cov_flag);
  
  dboot_t phys_res=cont_chir_linear_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,pars.adep_ml);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars]
		(double x,size_t ib)
		{return cont_chir_linear_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave());},
		bind(cont_chir_linear_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,pars.adep_ml),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));},
		ml_phys,phys_res,yaxis_label);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic quadratic fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_quad_ansatz(const Tpars &C0_ml,const Tpars &C1_ml,const Tpars &C2_ml,const Tm &ml,const Ta &a,const Tpars &adep,const Tpars &adep_ml)
{return C0_ml+C1_ml*ml+C2_ml*ml*ml+a*a*(adep+adep_ml*ml);}

//! perform the fit to the continuum limit
dboot_t cont_chir_quad_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,bool fix_adep_ml,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0,1},{0,1},{0,1},boot_fit);
  pars.add_adep_pars({0,1},{0,0},boot_fit);
  if(fix_adep_ml) boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_quad_ansatz(p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ac,p[pars.iadep],p[pars.iadep_ml]);},cov_flag);
  
  dboot_t phys_res=cont_chir_quad_ansatz(pars.C,pars.KPi,pars.KK,ml_phys,a_cont,pars.adep,pars.adep_ml);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars]
		(double x,size_t ib)
		{return cont_chir_quad_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave());},
		bind(cont_chir_quad_ansatz<dboot_t,double,double>,pars.C,pars.KPi,pars.KK,_1,a_cont,pars.adep,pars.adep_ml),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow));},
		ml_phys,phys_res,yaxis_label);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_Pi0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_Pi0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &Kpi0,
				   const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=(3+16*Cf04)*M2/den*log(M2);
  else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L)*pow(M2,0.5);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;

  return M2*Kpi0/(4*Cf04*den)*(1+chir_log/(4*Cf04)+K*M2/(4*Cf04*den))+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{2,10},{0.0,10},{80,10},boot_fit);
  pars.add_LEC_pars({5e-5,1e-5},{6,1},{0.0,10},boot_fit);
  //boot_fit.fix_par(pars.iadep_ml);
  //boot_fit.fix_par(pars.iL3dep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[chir_an](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon_Pi0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_Pi0(f0,B0,pars.C,pars.KPi,pars.KK,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"epsilon_Pi0: "<<phys_res.ave_err()<<", exp: 0.07"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_Pi0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_epsilon_Pi0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,_1,a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L)*pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\Pi0");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_K0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_K0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &Kk0,
				  const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    // M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    den=sqr(Tpars(4*M_PI*f0));
  
  // Tpars chir_log;
  // if(chir_an) chir_log=(3+16*Cf04)*M2Pi/den*log(M2Pi);
  // else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L)*pow(M2K,0.5);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;

  return M2K*Kk0/(4*Cf04*den)+disc_eff+fitted_FSE;
  // (1+chir_log/(4*Cf04)+K*M2Pi/(4*Cf04*den))
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{0.6,2},{0.0,10},{0.0,1},boot_fit);
  pars.add_LEC_pars({5e-5,1e-5},{6,1},{0.0,10},boot_fit);
  //boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[chir_an](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon_K0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_K0(f0,B0,pars.C,pars.KPi,pars.KK,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_K0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_epsilon_K0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L)*pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\K0");
  
  return phys_res;
}
