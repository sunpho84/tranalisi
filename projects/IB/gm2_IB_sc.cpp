#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr

template <class T,class Tm> T M2_fun(const T &B0,const Tm &aml)
{return 2*B0*aml;}

template <class T,class Tm> T M_fun(const T &B0,const Tm &aml)
{return sqrt(M2_fun(B0,aml));}

template <class T,class Tm> T xi_fun(const T &B0,const Tm &aml,const T &f0)
{return M2_fun(B0,aml)/sqr(T(4*M_PI*f0));}

template <class Tpars> Tpars FSE_LO(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
{return C*L3dep*xi*exp(-ML)/(ML);}

//! write an intestation line for asked ensemble
void write_ens_header(size_t iens)
{
  ens_data_t &ens=ens_data[iens];
  cout<<"----------------------- "<<iens<<" "<<ens.path<<" ---------------------"<<endl;
}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_LO(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
  Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  
  return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
}

//! perform the fit to the continuum limit of LO
dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
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
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_LO(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_LO(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,isyst]
		(double x,size_t ib)
		{return cont_chir_ansatz_LO<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
		bind(cont_chir_ansatz_LO<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{
		  auto aml=ext_data[idata].aml;
		  auto xi=xi_fun(B0,aml,f0),ML=M_fun(B0,aml)*ext_data[idata].L;
		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_LO(pars.C,pars.L3dep,xi,ML));},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={0.0,0.1};
  const ave_err_t K2Pi_guess={0,001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$a_\\mu";
  
  //set C
  ave_err_t C_guess;
  ave_err_t adep_guess;
  switch(im)
    {
    case istrange:
      C_guess=ave_err_t({5e-9,1e-9});
      adep_guess={0.0,0.1};
      break;
    case icharm:
      C_guess=ave_err_t({1.44e-9,1e-10});
      switch(case_of<c_cont>(isyst))
	{
	case 0:adep_guess={1.0,0.5};break;
	case 1:adep_guess={0,0.1};break;
	}
      break;
    default:
      CRASH("Unknwon mass %d",im);
    }
  
  cout<<"===================================================== LEADING ORDER ====================================================="<<endl;
  
  return cont_chir_fit_LO(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
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
dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
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
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,isyst]
		(double x,size_t ib)
		{return cont_chir_ansatz_QED<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
		bind(cont_chir_ansatz_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
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
      C_guess=ave_err_t({3e-12,0.5e-12});
      adep_guess={-1.5,0.6};
      KPi_guess={0.0,0.01};
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
dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
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
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_RAT(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_RAT(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,isyst]
		(double x,size_t ib)
		{return cont_chir_ansatz_RAT<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
		bind(cont_chir_ansatz_RAT<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_RAT(pars.C,pars.L3dep,ext_data[idata].L));},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={-0.001,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const string title="$$\\delta a_\\mu /a_\\mu";
  
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
  
  return cont_chir_fit_RAT(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

// template <class Tpars> Tpars FSE_m_rat(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
// {return C*L3dep*xi*exp(-ML)/(ML);}

// //! ansatz fit
// template <class Tpars,class Tm,class Ta>
// Tpars cont_chir_ansatz_m_rat(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t isyst)
// {
//   Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  
//   return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
// }

// //! perform the fit to the continuum limit of m_rat
// dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
// {
//   //set_printlevel(3);
  
//   boot_fit_t boot_fit;
//   size_t nbeta=a.size();
//   cont_chir_fit_pars_t pars(nbeta);
  
//   pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
//   pars.iC=boot_fit.add_fit_par(pars.C,"C_m_rat",C_guess.ave,C_guess.err);
//   pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
//   pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
//   pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
//   if(not FSE_an(isyst)) boot_fit.fix_par_to(pars.iL3dep,0.0);
//   boot_fit.fix_par_to(pars.iK2Pi,0.0);
//   if(chir_an(isyst)) boot_fit.fix_par_to(pars.iKPi,0.0);
//   if(cont_an(isyst)) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
//   cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
// 			 {return cont_chir_ansatz_m_rat(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
//   double a_cont=0;
//   dboot_t phys_res=cont_chir_ansatz_m_rat(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
//   cout<<"result: "<<phys_res.ave_err()<<endl;
  
//   plot_chir_fit(path,ext_data,pars,
// 		[&pars,isyst]
// 		(double x,size_t ib)
// 		{return cont_chir_ansatz_m_rat<double,double,double>
// 		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
// 		bind(cont_chir_ansatz_m_rat<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
// 		[&ext_data,&pars,&B0,&f0]
// 		(size_t idata,bool without_with_fse,size_t ib)
// 		{
// 		  auto aml=ext_data[idata].aml;
// 		  auto xi=xi_fun(B0,aml,f0),ML=M_fun(B0,aml)*ext_data[idata].L;
// 		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_m_rat(pars.C,pars.L3dep,xi,ML));
// 		},
// 		ml_phys,phys_res,title,beta_list);
  
//   return phys_res;
// }

// dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_ml_guess={0,0.001};
//   const ave_err_t KPi_guess={0.000,0.001};
//   const ave_err_t K2Pi_guess={-0.000,0.001};
//   const ave_err_t L3dep_guess={0,0.001};
//   const string title="$$m_{rat}-1";
  
//   ave_err_t C_guess;
//   ave_err_t adep_guess;
//   switch(im)
//     {
//     case istrange:
//       C_guess=ave_err_t({0.55,0.1});
//       adep_guess=ave_err_t({0.0,0.001});
//       break;
//     case icharm:
//       C_guess=ave_err_t({0.037,0.005});
//       adep_guess=ave_err_t({3.0,0.3});
//       break;
//     default:
//       CRASH("Unknwon mass %d",im);
//     }
  
//   cout<<"===================================================== MASS RATIO ====================================================="<<endl;
  
//   return cont_chir_fit_m_rat(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////////////

// template <class Tpars> Tpars FSE_dM(const Tpars &C,const Tpars &L3dep,double L)
// {return C*L3dep/L/L/L/L;}

// //! ansatz fit
// template <class Tpars,class Tm,class Ta>
// Tpars cont_chir_ansatz_dM(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t isyst)
// {
//   Tpars xi=xi_fun(B0,aml,f0);
  
//   return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_dM(C,L3dep,L);
// }

// //! perform the fit to the continuum limit of mass correction
// dboot_t cont_chir_fit_dM(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &title)
// {
//   //set_printlevel(3);
  
//   boot_fit_t boot_fit;
//   size_t nbeta=a.size();
//   cont_chir_fit_pars_t pars(nbeta);
  
//   pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
//   pars.iC=boot_fit.add_fit_par(pars.C,"C_dM",C_guess.ave,C_guess.err);
//   pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
//   pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
//   pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
//   if(not FSE_an(isyst)) boot_fit.fix_par_to(pars.iL3dep,0.0);
//   boot_fit.fix_par_to(pars.iK2Pi,0.0);
//   if(chir_an(isyst)) boot_fit.fix_par_to(pars.iKPi,0.0);
//   if(cont_an(isyst)) boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
//   cont_chir_fit_minimize(ext_data,pars,boot_fit,1.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
// 			 {return cont_chir_ansatz_dM(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
//   double a_cont=0;
//   dboot_t phys_res=cont_chir_ansatz_dM(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
//   cout<<"result: "<<phys_res.ave_err()<<endl;
  
//   plot_chir_fit(path,ext_data,pars,
// 		[&pars,isyst]
// 		(double x,size_t ib)
// 		{return cont_chir_ansatz_dM<double,double,double>
// 		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
// 		bind(cont_chir_ansatz_dM<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
// 		[&ext_data,&pars,&B0,&f0]
// 		(size_t idata,bool without_with_fse,size_t ib)
// 		{return dboot_t(ext_data[idata].wfse/pars.fit_a[ib]-without_with_fse*FSE_dM(pars.C,pars.L3dep,ext_data[idata].L));},
// 		ml_phys,phys_res,title,beta_list);
  
//   return phys_res;
// }

// dboot_t cont_chir_fit_dM(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_ml_guess={0,0.001};
//   const ave_err_t KPi_guess={0.000,0.001};
//   const ave_err_t K2Pi_guess={-0.000,0.001};
//   const ave_err_t L3dep_guess={0,0.001};
//   const string title="dm";
  
//   ave_err_t C_guess;
//   ave_err_t adep_guess;
//   switch(im)
//     {
//     case istrange:
//       C_guess=ave_err_t({0.55,0.1});
//       adep_guess=ave_err_t({0.0,0.001});
//       break;
//     case icharm:
//       C_guess=ave_err_t({0.037,0.005});
//       adep_guess=ave_err_t({3.0,0.3});
//       break;
//     default:
//       CRASH("Unknwon mass %d",im);
//     }
  
//   cout<<"===================================================== MASS CORRECTION ====================================================="<<endl;
  
//   return cont_chir_fit_dM(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,title);
// }

///////////////////////////////////////////////////////////////////////////////////////////////

void perform_analysis(const dbvec_t &data,const string &name)
{
  double ave=0,err=0,wave=0,werr=0,weight=0;
  for(auto &d : data)
    {
      double a=d.ave(),e=d.err(),w=1/sqr(e);
      ave+=a;
      err+=e;
      wave+=a*w;
      werr+=e*w;
      weight+=w;
    }
  ave/=data.size();
  err/=data.size();
  wave/=weight;
  werr/=weight;
  vector<double> syst=syst_analysis_sep_bis(data.ave_err(),ind_syst);
  double tot=0;
  for(auto &s : syst) tot+=s*s;
  tot=sqrt(tot);
  
  cout<<" ================== "<<name<<" ================== "<<endl;
  cout<<" Aver:\t"<<ave<<endl;
  cout<<"  Weig:\t"<<wave<<endl;
  vector<pair<basic_string<char>,double>> res({{"Stat",err},{" Weig",werr}});
  for(size_t isyst=0;isyst<ind_syst.rank();isyst++) res.push_back({ind_syst.name(isyst),syst[isyst]});
  res.push_back({"TotSy",tot});
  res.push_back({"Tot",sqrt(sqr(tot)+sqr(err))});
  for(auto &el : res)
    {
      streamsize prec=cout.precision();
      cout<<" "<<el.first.substr(0,6)<<":\t"<<el.second<<"\t=\t";
      cout.precision(3);
      cout<<fabs(el.second/ave)*100<<"% "<<endl;
      tot+=sqr(el.second);
      cout.precision(prec);
    }
  
  // for(auto x : data.ave_err()) cout<<" "<<x<<endl;
}

int main(int narg,char **arg)
{
  // vector<ave_err_t> v({{-0.001518092727373818,0},{-0.001563828393244764,0},{-0.001570134136155057,0},{-0.001555232660649148,0},{-0.001385082821679269,0},{-0.001473247256814399,0},{-0.001367408317152487,0},{-0.001448039412439519,0}});
  
  // vector<ave_err_t> v1({
  // 			{-0.001518092727373818,0},
  // 			{-0.001518092727373818,0},
  // 			{-0.001518092727373818,0},
  // 			{-0.001518092727373818,0},
  // 			{-0.001385082821679269,0},
  // 			{-0.001385082821679269,0},
  // 			{-0.001385082821679269,0},
  // 			{-0.001385082821679269,0}}
  //   );
  // syst_analysis_sep(v1);
   
  //  vector<double> res=syst_analysis_sep_bis(v1,{2,2,2});
  //  cout<<endl;
   
  //  for(auto x : res) cout<<"     "<<x<<endl;
  //  cout<<endl;
  
  gm2_initialize(narg,arg);
  
  vector<djvec_t> jPP_LO(nens_used),jPP_QED(nens_used),jVV_LO(nens_used),jVV_QED(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      djack_t deltam_cr=compute_deltam_cr(ens,ilight);
      
      //load LO for PP
      jPP_LO[iens]=read_PP("00",ens,im,1,RE);
      jPP_QED[iens]=read_QED_PP(ens,1,im,deltam_cr,jPP_LO[iens]);
      
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
  
  //fit 2pts
  index_t ind_2pts_fit({{"NFitRanges",nfit_range_variations},{"Ens",nens_used}});
  size_t n2pts_fit_max=ind_2pts_fit.max();
  djvec_t jZ_V(n2pts_fit_max),jM_V(n2pts_fit_max),jA_V(n2pts_fit_max),jSL_V(n2pts_fit_max);
  djvec_t jZ_P(n2pts_fit_max),jM_P(n2pts_fit_max),jA_P(n2pts_fit_max),jSL_P(n2pts_fit_max);
  
  //read if avail
  for(auto &obj : {&jZ_P,&jM_P,&jA_P,&jSL_P,&jZ_V,&jM_V,&jA_V,&jSL_V})
    if(an_mode!=compute_everything) obj->bin_read(results_in);
  
  cout<<" ********************************************* fitting 2pts *******************************************"<<endl;
  for(size_t ind=0;ind<n2pts_fit_max;ind++)
    {
      vector<size_t> comp=ind_2pts_fit(ind);
      size_t ifit_range=comp[0],iens=comp[1];
      write_ens_header(iens);
      cout<<"Flags: ifit_range="<<ifit_range<<", iens="<<iens<<endl;
      
      //set malus for fitting
      int malus_fitting=0;
      if(ifit_range) malus_fitting=1;
      
      ens_data_t &ens=ens_data[iens];
      size_t TH=ens.T/2;
      string ens_qpath=ens.path+"/plots_"+qname[im];
      
      if(an_mode==compute_everything)
	{
	  two_pts_with_ins_ratio_fit(jZ_V[ind],jM_V[ind],jA_V[ind],jSL_V[ind],jVV_LO[iens],jVV_QED[iens],TH,ens.tmin[im]+malus_fitting,ens.tmax[im],
				     combine("%s/VV_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/VV_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	  two_pts_with_ins_ratio_fit(jZ_P[ind],jM_P[ind],jA_P[ind],jSL_P[ind],jPP_LO[iens],jPP_QED[iens],TH,ens.tmin[im]+malus_fitting,ens.tmax[im],
				     combine("%s/PP_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/PP_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	}
      
      cout<<"M_V: "<<jM_V[ind].ave_err()<<endl;
      cout<<"M_P: "<<jM_P[ind].ave_err()<<endl;
    }
  
  //write if computed
  for(auto &obj : {&jZ_P,&jM_P,&jA_P,&jSL_P,&jZ_V,&jM_V,&jA_V,&jSL_V})
    if(an_mode==compute_everything) obj->bin_write(results_out);
  
  //perform the integrations
  cout<<endl;
  cout<<" ********************************************* integrating *******************************************"<<endl;
  index_t ind_integr({{"Input",ninput_an},{"Cont",ncont_extrap},{"FitRange",nfit_range_variations},{"Ens",nens_used}});
  size_t nintegr_max=ind_integr.max();
  dbvec_t resc_a(nintegr_max);
  dbvec_t  LO_correl(nintegr_max), LO_remaind(nintegr_max), LO(nintegr_max);
  dbvec_t QED_correl(nintegr_max),QED_remaind(nintegr_max),QED(nintegr_max);
  
  //read if avail
  for(auto &obj : {&resc_a,&LO_correl,&LO_remaind,&QED_correl,&QED_remaind})
    if(an_mode!=compute_everything) obj->bin_read(results_in);
  
  for(size_t ind=0;ind<nintegr_max;ind++)
    {
      vector<size_t> comp=ind_integr(ind);
      size_t input_an_id=comp[0],icont_extrap=comp[1],ifit_range=comp[2],iens=comp[3];
      ens_data_t &ens=ens_data[iens];
      string ens_qpath=ens.path+"/plots_"+qname[im];
      size_t ib=ens.ib,TH=ens.T/2;
      write_ens_header(iens);
      cout<<"Flags: input_an_id="<<input_an_id<<", icont_extrap="<<icont_extrap<<", ifit_range="<<ifit_range<<", iens="<<iens<<endl;
      
      size_t i2pts=ind_2pts_fit({ifit_range,iens});
      
      //set bootstrap
      bi=jack_index[input_an_id][ens.iult];
      prepare_az(input_an_id);
      dbvec_t VV_LO(bi,jVV_LO[iens]),PP_LO(bi,jPP_LO[iens]);
      dbvec_t VV_QED(bi,jVV_QED[iens]),PP_QED(bi,jPP_QED[iens]);
      dboot_t Z_V(bi,jZ_V[i2pts]),M_V(bi,jM_V[i2pts]),A_V(bi,jA_V[i2pts]),SL_V(bi,jSL_V[i2pts]);
      
      if(an_mode==compute_everything)
	{
	  //lattice spacing obtained from V
	  switch(icont_extrap)
	    {
	    case 0:resc_a[ind]=1/lat_par[input_an_id].ainv[ib];break;
	    case 1:resc_a[ind]=M_V/M_V_phys[im];break;
	    }
	  cout<<"a: "<<resc_a[ind].ave_err()<<endl;
	  
	  size_t upto=TH-DT;
	  LO_correl[ind]=integrate_corr_times_kern_up_to(VV_LO,ens.T,resc_a[ind],im,upto)*sqr(Za[ib]);
	  LO_remaind[ind]=integrate_LO_reco_from(Z_V,M_V,resc_a[ind],im,upto)*sqr(Za[ib]);
	  
	  QED_correl[ind]=integrate_corr_times_kern_up_to(VV_QED,ens.T,resc_a[ind],im,upto)*sqr(Za[ib]);
	  QED_remaind[ind]=integrate_QED_reco_from(A_V,Z_V,SL_V,M_V,resc_a[ind],im,upto)*sqr(Za[ib]);
	  
	  index_t ind_rest({{"Input",ninput_an},{"Cont",ncont_extrap},{"FitRange",nfit_range_variations}});
	  size_t irest=ind_rest({input_an_id,icont_extrap,ifit_range});
	  
	  compare_LO_num_reco(combine("%s/kern_LO_num_reco_%zu.xmg",ens_qpath.c_str(),irest),VV_LO,Z_V,M_V,resc_a[iens]);
	  compare_QED_num_reco(combine("%s/kern_QED_num_reco_%zu.xmg",ens_qpath.c_str(),irest),VV_QED,A_V,Z_V,SL_V,M_V,resc_a[iens]);
	}
      
      LO[ind]=LO_correl[ind]+LO_remaind[ind];
      QED[ind]=QED_correl[ind]+QED_remaind[ind];
      
      cout<<"amu: "<< LO_correl[ind].ave_err()<<" + "<<LO_remaind[ind].ave_err()<<" = "<<LO[ind].ave_err()<<endl;
      cout<<"amu_QED: "<<QED_correl[ind].ave_err()<<" + "<<QED_remaind[ind].ave_err()<<" = "<<QED[ind].ave_err()<<endl;
    }
  
  //write if computed
  for(auto &obj : {&resc_a,&LO_correl,&LO_remaind,&QED_correl,&QED_remaind})
    if(an_mode==compute_everything) obj->bin_write(results_out);
  
  //loop over analysis flags and input scale determination
  dbvec_t cLO(ind_syst.max());
  dbvec_t cQED(ind_syst.max());
  dbvec_t cRAT(ind_syst.max());
  //dbvec_t c_m_rat,c_dMP,c_dMV;
  for(size_t isyst=0;isyst<ind_syst.max();isyst++)
    {
      vector<size_t> comp=ind_syst(isyst);
      cout<<"Analysis: "<<isyst<<endl;
      size_t input_an_id=case_of<c_input>(isyst);
      size_t chir_an_id=case_of<c_chir>(isyst);
      size_t FSE_an_id=case_of<c_FSE>(isyst);
      size_t cont_an_id=case_of<c_cont>(isyst);
      size_t fit_range_an_id=case_of<c_fit_range>(isyst);
      cout<<"InputAn:\t"<<input_an_id<<endl;
      cout<<"ChirAn: \t"<<chir_an_id<<endl;
      cout<<"FSE:    \t"<<FSE_an_id<<endl;
      cout<<"ContAn: \t"<<cont_an_id<<endl;
      cout<<"FitRange:\t"<<fit_range_an_id<<endl;
      
      prepare_az(input_an_id);
      dboot_t &f0=lat_par[input_an_id].f0;
      dboot_t &B0=lat_par[input_an_id].B0;
      
      if(case_of<c_FSE>(isyst))
	{
	  grace::default_color_scheme={grace::RED,grace::RED,grace::RED, grace::BLUE,grace::BLUE, grace::GREEN4,grace::VIOLET};
	  grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
	}
      else
	{
	  grace::default_color_scheme={grace::RED, grace::BLUE, grace::GREEN4,grace::VIOLET};
	  grace::default_symbol_scheme={grace::DIAMOND,grace::DIAMOND,grace::DIAMOND};
	}
      
      vector<cont_chir_fit_data_t> data_LO;
      vector<cont_chir_fit_data_t> data_QED;
      vector<cont_chir_fit_data_t> data_RAT;
      // vector<cont_chir_fit_data_t> data_dMP;
      // vector<cont_chir_fit_data_t> data_dMV;
      // vector<cont_chir_fit_data_t> data_m_rat;
      for(size_t iens=0;iens<ens_data.size();iens++)
	if(case_of<c_FSE>(isyst) or ens_data[iens].use_for_L)
	  {
	    write_ens_header(iens);
	    ens_data_t &ens=ens_data[iens];
	    string ens_qpath=ens.path+"/plots_"+qname[im];
	    
	    size_t iintegr=ind_integr({input_an_id,cont_an_id,fit_range_an_id,iens});
	    
	    double rat_perturb=-0.01835*sqr(eq[im]);
	    dboot_t RAT=QED[iintegr]/LO[iintegr]+rat_perturb;
	    cout<<" Ratio: "<<RAT.ave_err()<<endl;
	    
	    //passing dummy aux mass
	    dboot_t dum;
	    dum=0.0;
	    data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,LO[iintegr],LO[iintegr]));
	    data_QED.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,QED[iintegr],QED[iintegr]));
	    data_RAT.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,RAT,RAT));
	    
	    // //cc ratio
	    // dboot_t m_rat=M_V/M_P[iens]-1.0;
	    // data_m_rat.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,m_rat,m_rat));
	    
	    // //mass corrections
	    // data_dMP.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,SL_P,SL_P));
	    // data_dMV.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,SL_V,SL_V));
	  }
      
      cLO[isyst]=cont_chir_fit_LO(alist,zlist,f0,B0,data_LO,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_LO_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      cQED[isyst]=cont_chir_fit_QED(alist,zlist,f0,B0,data_QED,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_QED_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      cRAT[isyst]=cont_chir_fit_RAT(alist,zlist,f0,B0,data_RAT,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_RAT_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      
      // c_m_rat[iai]=cont_chir_fit_m_rat(alist,zlist,f0,B0,data_m_rat,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_m_rat_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
      
      // c_dMP[iai]=cont_chir_fit_dM(alist,zlist,f0,B0,data_dMP,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_dMP_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
      // c_dMV[iai]=cont_chir_fit_dM(alist,zlist,f0,B0,data_dMV,lat_par[input_an_id].ml,combine("plots_%s/cont_chir_dMV_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
    }
  
  perform_analysis(cLO,"LO");
  perform_analysis(cQED,"QED");
  perform_analysis(cRAT,"RAT");
  
  // perform_analysis(c_m_rat,"m_rat");
  
  return 0;
}
