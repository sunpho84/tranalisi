#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr
const int im=istrange;
//const int im=icharm;

template <class T,class Tm> T M2_fun(const T &B0,const Tm &aml)
{return 2*B0*aml;}

template <class T,class Tm> T M_fun(const T &B0,const Tm &aml)
{return sqrt(M2_fun(B0,aml));}

template <class T,class Tm> T xi_fun(const T &B0,const Tm &aml,const T &f0)
{return M2_fun(B0,aml)/sqr(T(4*M_PI*f0));}

template <class Tpars> Tpars FSE_LO(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
{return C*L3dep*xi*exp(-ML)/(ML);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_LO(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  
  return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
}

//! perform the fit to the continuum limit of LO
dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
{
  //set_printlevel(3);
  //if(an_flag==1) distr_fit_debug=1;
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_LO",C_guess.ave,C_guess.err);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
  if(not FSE_an(an_flag)) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)) boot_fit.fix_par_to(pars.iKPi,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_LO(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_LO(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_LO<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_LO<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{
		  auto aml=ext_data[idata].aml;
		  auto xi=xi_fun(B0,aml,f0),ML=M_fun(B0,aml)*ext_data[idata].L;
		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_LO(pars.C,pars.L3dep,xi,ML));},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={0.0,0.1};
  const ave_err_t K2Pi_guess={0,001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$a_\\mu";
  
  //set C
  ave_err_t C_guess;
  switch(im)
    {
    case istrange:
      C_guess=ave_err_t({5e-9,1e-9});
      break;
    case icharm:
      C_guess=ave_err_t({1.44e-9,1e-10});
      break;
    default:
      CRASH("Unknwon mass %d",im);
    }
  
  //set adep
  ave_err_t adep_guess;
  if(cont_an(an_flag)) adep_guess={0.5,0.08};
  else                 adep_guess={1,0.5};
  
  cout<<"===================================================== LEADING ORDER ====================================================="<<endl;
  
  return cont_chir_fit_LO(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tpars> Tpars FSE_QED(const Tpars &C,const Tpars &L3dep,const double L)
{return C*L3dep/L/L/L;}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars xi=xi_fun(B0,aml,f0);
  return C*(1.0+Kpi*xi+Kpi*xi*xi+a*a*adep)+FSE_QED(C,L3dep,L);
}

//! perform the fit to the continuum limit of QED
dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_QED",C_guess.ave,C_guess.err);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
  if(not FSE_an(an_flag)) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_QED<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_QED(pars.C,pars.L3dep,ext_data[idata].L));},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$a^{QED}_\\mu";
  
  ave_err_t C_guess;
  ave_err_t adep_guess;
  ave_err_t KPi_guess;
  switch(im)
    {
    case istrange:
      C_guess=ave_err_t({1e-10,1e-12});
      adep_guess={-1.37,0.8};
      KPi_guess={0.0,0.1};
      break;
    case icharm:
      C_guess=ave_err_t({-7e-12,2e-12});
      adep_guess={-6,0.2};
      KPi_guess={-0.001,0.001};
      break;
    default:
      CRASH("Unknwon mass %d",im);
    }
  
  cout<<"===================================================== QED CORR ====================================================="<<endl;
  
  return cont_chir_fit_QED(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tpars> Tpars FSE_RAT(const Tpars &C,const Tpars &L3dep,const double L)
{return C*L3dep/L/L/L;}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_RAT(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars xi=xi_fun(B0,aml,f0);
  return C*(1.0+Kpi*xi+Kpi*xi*xi+a*a*adep)+FSE_RAT(C,L3dep,L);
}

//! perform the fit to the continuum limit of RAT
dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_RAT",C_guess.ave,C_guess.err);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
  if(not FSE_an(an_flag)) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_RAT(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_RAT(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_RAT<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_RAT<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_RAT(pars.C,pars.L3dep,ext_data[idata].L));},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={-0.001,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$\\delta a_\\mu/a_\\mu";
  
  ave_err_t C_guess;
  ave_err_t adep_guess;
  switch(im)
    {
    case istrange:
      C_guess=ave_err_t({-0.0018,0.0004});
      adep_guess={-1,0.2};
  break;
    case icharm:
      C_guess=ave_err_t({-0.013,0.002});
      adep_guess={-2,0.3};
      break;
    default:
      CRASH("Unknwon mass %d",im);
    }
  
  cout<<"===================================================== RATIO ====================================================="<<endl;
  
  return cont_chir_fit_RAT(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tpars> Tpars FSE_m_rat(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
{return C*L3dep*xi*exp(-ML)/(ML);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_m_rat(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  
  return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
}

//! perform the fit to the continuum limit of m_rat
dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_m_rat",C_guess.ave,C_guess.err);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
  if(not FSE_an(an_flag)) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(chir_an(an_flag)) boot_fit.fix_par_to(pars.iKPi,0.0);
  if(cont_an(an_flag)) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_m_rat(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_m_rat(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_m_rat<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_m_rat<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{
		  auto aml=ext_data[idata].aml;
		  auto xi=xi_fun(B0,aml,f0),ML=M_fun(B0,aml)*ext_data[idata].L;
		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_m_rat(pars.C,pars.L3dep,xi,ML));
		},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={0.000,0.001};
  const ave_err_t K2Pi_guess={-0.000,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$m_{rat}-1";
  
  ave_err_t C_guess;
  ave_err_t adep_guess;
  switch(im)
    {
    case istrange:
      C_guess=ave_err_t({0.55,0.1});
      adep_guess=ave_err_t({0.0,0.001});
      break;
    case icharm:
      C_guess=ave_err_t({0.037,0.005});
      adep_guess=ave_err_t({3.0,0.3});
      break;
    default:
      CRASH("Unknwon mass %d",im);
    }
  
  cout<<"===================================================== MASS RATIO ====================================================="<<endl;
  
  return cont_chir_fit_m_rat(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
}

///////////////////////////////////////////////////////////////////////////////////////////////

void perform_analysis(const dbvec_t &data,const string &name)
{
  vector<ave_err_t> ave=ave_analyses(data);
  cout<<" ================== "<<name<<" ================== "<<endl;
  cout<<stat_analysis(ave)<<" "<<syst_analysis(ave)<<endl;
   for(size_t i=0;i<ave.size();i++)
     cout<<"  "<<ave[i]<<endl;
   syst_analysis_sep(ave);
}

int main(int narg,char **arg)
{
  gm2_initialize(narg,arg);
  
  vector<djvec_t> jPP_LO(nens_used),jVV_LO(nens_used),jVV_QED(nens_used);
  djvec_t M_P(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      size_t TH=ens.T/2;
      djack_t deltam_cr=compute_deltam_cr(ens,ilight);
      
      //load LO for PP
      jPP_LO[iens]=read_PP("00",ens,im,1,RE);
      djack_t Z_P;
      two_pts_fit(Z_P,M_P[iens],jPP_LO[iens],TH,ens.tmin[im],ens.tmax[im],combine("%s/plots/PP_LO_emass.xmg",ens.path.c_str()));
      
      //load LO and QED for VV
      jVV_LO[iens]=read_VV("00",ens,im,1,RE);
      jVV_QED[iens]=read_QED_VV(ens,1,im,deltam_cr,jVV_LO[iens]);
    }
  
  vector<string> beta_list={"1.90","1.95","2.10"};
  
  string res_path=combine("integr_results_%s",qname[im].c_str());
  enum an_mode_t{compute_everything,read_integrals};
  an_mode_t an_mode=compute_everything;
  
  //files for in/out
  raw_file_t results_in;
  raw_file_t results_out;
  
  //check if the results exists
  if(file_exists(res_path))
    {
      an_mode=read_integrals;
      results_in.open(res_path,"r");
      cout<<"Reading integrals from the file "<<res_path<<endl;
    }
  else
    {
      an_mode=compute_everything;
      results_out.open(res_path,"w");
      cout<<"File "<<res_path<<" absent, computing everything"<<endl;
    }
  
  //loop over analysis flags and input scale determination
  dbvec_t cLO(ninput_an*nan_syst);
  dbvec_t cQED(ninput_an*nan_syst);
  dbvec_t cRAT(ninput_an*nan_syst);
  dbvec_t c_m_rat(ninput_an*nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	prepare_az(input_an_id);
	dboot_t &f0=lat_par[input_an_id].f0;
	dboot_t &B0=lat_par[input_an_id].B0;
	
	cout<<" ANALYSIS flag: "<<an_flag<<" , input_an_id: "<<input_an_id<<endl;
	
	if(FSE_an(an_flag))
	  {
	    grace::default_color_scheme={grace::RED,grace::RED,grace::RED, grace::BLUE,grace::BLUE, grace::GREEN4,grace::VIOLET};
	    grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
	  }
	else
	  {
	    grace::default_color_scheme={grace::RED, grace::BLUE, grace::GREEN4,grace::VIOLET};
	    grace::default_symbol_scheme={grace::DIAMOND,grace::DIAMOND,grace::DIAMOND};
	  }
	
	dboot_t LO_correl,LO_remaind,LO,QED_correl,QED_remaind,QED;
	
	vector<cont_chir_fit_data_t> data_LO;
	vector<cont_chir_fit_data_t> data_QED;
	vector<cont_chir_fit_data_t> data_RAT;
	vector<cont_chir_fit_data_t> data_m_rat;
	for(size_t iens=0;iens<ens_data.size();iens++)
	  if(FSE_an(an_flag) or ens_data[iens].use_for_L)
	    {
	      ens_data_t &ens=ens_data[iens];
	      cout<<"----------------------- "<<iens<<" "<<ens.path<<" ---------------------"<<endl;
	      int ib=ens.ib;
	      size_t TH=ens.T/2;
	      
	      //set bootstrap
	      bi=jack_index[input_an_id][ens.iult];
	      dbvec_t VV_LO(bi,jVV_LO[iens]);
	      dbvec_t VV_QED(bi,jVV_QED[iens]);
	      
	      dboot_t Z_V,M_V,A_V,SL_V;
	      if(an_mode==compute_everything)
		two_pts_with_ins_ratio_fit(Z_V,M_V,A_V,SL_V,VV_LO,VV_QED,TH,ens.tmin[im],ens.tmax[im],
					   combine("%s/plots/VV_LO.xmg",ens.path.c_str()),combine("%s/plots/VV_QED.xmg",ens.path.c_str()));
	      //read or write
	      for(auto &obj : {&Z_V,&M_V,&A_V,&SL_V})
	       	if(an_mode==compute_everything) obj->bin_write(results_out);
		else 	                        obj->bin_read(results_in);
	      
	      //lattice spacing obtained from V
	      dboot_t resc_a;
	      if(cont_an(an_flag)) resc_a=M_V/M_V_phys[im];
	      else                 resc_a=1/lat_par[input_an_id].ainv[ib];
	      
	      if(an_mode==compute_everything)
		{
#pragma omp parallel sections
		  {
#pragma omp section
		    {
		      size_t upto=TH-DT;
		      LO_correl=integrate_corr_times_kern_up_to(VV_LO,ens.T,resc_a,im,upto)*sqr(Za[ib]);
		      LO_remaind=integrate_LO_reco_from(Z_V,M_V,resc_a,im,upto)*sqr(Za[ib]);
		    }
#pragma omp section
		    {
		      size_t upto=TH-DT;
		      QED_correl=integrate_corr_times_kern_up_to(VV_QED,ens.T,resc_a,im,upto)*sqr(Za[ib]);
		      QED_remaind=integrate_QED_reco_from(A_V,Z_V,SL_V,M_V,resc_a,im,upto)*sqr(Za[ib]);
		    }
		  }
		  
		  LO_correl.bin_write(results_out);
		  LO_remaind.bin_write(results_out);
		  QED_correl.bin_write(results_out);
		  QED_remaind.bin_write(results_out);
		}
	      else
		{
		  LO_correl.bin_read(results_in);
		  LO_remaind.bin_read(results_in);
		  QED_correl.bin_read(results_in);
		  QED_remaind.bin_read(results_in);
		}
	      
	      LO=LO_correl+LO_remaind;
	      QED=QED_correl+QED_remaind;
	      
	      compare_LO_num_reco(combine("%s/plots/kern_LO_num_reco.xmg",ens.path.c_str()),VV_LO,Z_V,M_V,resc_a);
	      cout<<"amu: "<< LO_correl.ave_err()<<" + "<<LO_remaind.ave_err()<<" = "<<LO.ave_err()<<endl;
	      
	      compare_QED_num_reco(combine("%s/plots/kern_QED_num_reco.xmg",ens.path.c_str()),VV_QED,A_V,Z_V,SL_V,M_V,resc_a);
	      cout<<"amu_QED: "<<QED_correl.ave_err()<<" + "<<QED_remaind.ave_err()<<" = "<<QED.ave_err()<<endl;
	      
	      double rat_perturb=-0.01835*sqr(eq[im]);
	      dboot_t RAT=QED/LO+rat_perturb;
	      cout<<" Ratio: "<<RAT.ave_err()<<endl;
	      
	      //to avoid D
	      dboot_t dum;
	      dum=0.0;
	      data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,LO,LO));
	      data_QED.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,QED,QED));
	      data_RAT.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,RAT,RAT));
	      
	      //cc ratio
	      dboot_t m_rat=M_V/dboot_t(bi,M_P[iens])-1.0;
	      data_m_rat.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,m_rat,m_rat));
	    }
	
	int iai=ind_an({input_an_id,an_flag});
	
	cLO[iai]=cont_chir_fit_LO(alist,zlist,f0,B0,data_LO,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_LO_flag%zu_an%zu.xmg",qname[im].c_str(),an_flag,input_an_id),an_flag,false,beta_list);
	cQED[iai]=cont_chir_fit_QED(alist,zlist,f0,B0,data_QED,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_QED_flag%zu_an%zu.xmg",qname[im].c_str(),an_flag,input_an_id),an_flag,false,beta_list);
	cRAT[iai]=cont_chir_fit_RAT(alist,zlist,f0,B0,data_RAT,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_RAT_flag%zu_an%zu.xmg",qname[im].c_str(),an_flag,input_an_id),an_flag,false,beta_list);
	
	c_m_rat[iai]=cont_chir_fit_m_rat(alist,zlist,f0,B0,data_m_rat,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_m_rat_flag%zu_an%zu.xmg",qname[im].c_str(),an_flag,input_an_id),an_flag,false,beta_list);
      }
  
   perform_analysis(cLO,"LO");
   perform_analysis(cQED,"QED");
   perform_analysis(cRAT,"RAT");
   
   perform_analysis(c_m_rat,"m_rat");
   
   return 0;
}
