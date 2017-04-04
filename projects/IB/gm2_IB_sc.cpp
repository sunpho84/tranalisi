#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr

enum an_mode_t{compute_everything,read_integrals};
an_mode_t an_mode=compute_everything;

template <class Tpars> Tpars FSE_LO(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
{return C*L3dep*xi*exp(-ML)/(ML);}

//return the tmin for the fit
size_t tmin_fit(size_t iens,size_t ifit_range)
{
  //set malus for fitting
  int malus_fitting=0;
  if(ifit_range) malus_fitting=2;
  
  ens_data_t &ens=ens_data[iens];
  
  return ens.tmin[im]+malus_fitting;
}

//! write an intestation line for asked ensemble
void write_ens_header(size_t iens)
{
  ens_data_t &ens=ens_data[iens];
  cout<<"----------------------- "<<iens<<" "<<ens.path<<" ---------------------"<<endl;
}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_LO(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t &isyst)
{
  Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
}

//! perform the fit to the continuum limit of LO
dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &yaxis_title)
{
  //set_printlevel(3);
  //if(isyst==1) distr_fit_debug=1;
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_LO",C_guess.ave(),C_guess.err());
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  
  if(case_of<c_FSE>(isyst)%2==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst]
			 (const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_LO(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],
						     p[pars.iadep_ml],L,p[pars.iL3dep],isyst);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_LO(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,isyst]
		(double x,size_t ib)
		{return cont_chir_ansatz_LO<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),
		     pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),isyst);},
		bind(cont_chir_ansatz_LO<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,
		     pars.adep,pars.adep_ml,inf_vol,pars.L3dep,isyst),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{
		  dboot_t aml=ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib];
		  dboot_t xi=xi_fun(pars.fit_B0,aml,pars.fit_f0),ML=M_fun(pars.fit_B0,aml)*ext_data[idata].L;
		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_LO(pars.C,pars.L3dep,xi,ML));},
		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  return phys_res;
}

dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={0.0,0.1};
  const ave_err_t K2Pi_guess={0,001};
  const ave_err_t L3dep_guess={0,0.001};
  const string yaxis_title="$$a_\\mu";
  
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
  return cont_chir_fit_LO(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,yaxis_title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tpars> Tpars FSE_QED(const Tpars &C,const Tpars &L3dep,const double &L,const Tpars &powL)
{return C*L3dep*pow(L,powL);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const Tpars &powL)
{
  Tpars xi=xi_fun(B0,aml,f0);
  return C*(1.0+Kpi*xi+Kpi*xi*xi+a*a*adep)+FSE_QED(C,L3dep,L,powL);
}

//! perform the fit to the continuum limit of QED
dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const ave_err_t &powL_guess,const string &yaxis_title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_QED",C_guess.ave(),C_guess.err());
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  pars.iL4dep=boot_fit.add_fit_par(pars.L4dep,"powL",powL_guess.ave(),powL_guess.err());
  boot_fit.fix_par_to(pars.iL4dep,-2.0);
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],p[pars.iL4dep]);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,pars.L4dep);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars]
		(double x,size_t ib)
		{return cont_chir_ansatz_QED<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),
		     inf_vol,pars.L3dep.ave(),pars.L4dep.ave());},
		bind(cont_chir_ansatz_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,
		     inf_vol,pars.L3dep,pars.L4dep),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_QED(pars.C,pars.L3dep,ext_data[idata].L,pars.L4dep));},
		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  return phys_res;
}

dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const ave_err_t powL_guess={-2,0.001};
  const string yaxis_title="$$a^{QED}_\\mu";
  
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
  
  return cont_chir_fit_QED(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,powL_guess,yaxis_title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class Tpars> Tpars FSE_RAT(const Tpars &C,const Tpars &L3dep,const double &L,const Tpars &powL)
{return FSE_QED(C,L3dep,L,powL);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_RAT(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const Tpars &powL)
{
  Tpars xi=xi_fun(B0,aml,f0);
  return C*(1.0+Kpi*xi+Kpi*xi*xi+a*a*adep)+FSE_RAT(C,L3dep,L,powL);
}

//! perform the fit to the continuum limit of RAT
dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const ave_err_t &powL_guess,const string &yaxis_title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_RAT",C_guess.ave(),C_guess.err());
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  pars.iL4dep=boot_fit.add_fit_par(pars.L4dep,"powL",powL_guess.ave(),powL_guess.err());
  boot_fit.fix_par_to(pars.iL4dep,-2.0);

  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst]
			 (const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_RAT(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,
						      p[pars.iL3dep],p[pars.iL4dep]);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_RAT(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,pars.L4dep);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars]
		(double x,size_t ib)
		{return cont_chir_ansatz_RAT<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),
		     inf_vol,pars.L3dep.ave(),pars.L4dep.ave());},
		bind(cont_chir_ansatz_RAT<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,
		     inf_vol,pars.L3dep,pars.L4dep),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_RAT(pars.C,pars.L3dep,ext_data[idata].L,pars.L4dep));},
		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  return phys_res;
}

dboot_t cont_chir_fit_RAT(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t KPi_guess={-0.001,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const ave_err_t powL_guess={-2,0.001};
  const string yaxis_title="$$\\delta a_\\mu /a_\\mu";
  
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
  
  return cont_chir_fit_RAT(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,powL_guess,yaxis_title);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

// template <class Tpars> Tpars FSE_m_rat(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &ML)
// {return C*L3dep*xi*exp(-ML)/(ML);}

// //! ansatz fit
// template <class Tpars,class Tm,class Ta>
// Tpars cont_chir_ansatz_m_rat(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t &isyst)
// {
//   Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml),ML=M*L;
  
//   return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,ML);
// }

// //! perform the fit to the continuum limit of m_rat
// dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &yaxis_title)
// {
//   //set_printlevel(3);
  
//   boot_fit_t boot_fit;
//   size_t nbeta=a.size();
//   cont_chir_fit_pars_t pars(nbeta);
  
//   pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
//   pars.iC=boot_fit.add_fit_par(pars.C,"C_m_rat",C_guess.ave(),C_guess.err());
//   pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
//   pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
//   pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  
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
//		  dboot_t aml=ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib];
// 		  dboot_t xi=xi_fun(B0,aml,f0),ML=M_fun(B0,aml)*ext_data[idata].L;
// 		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_m_rat(pars.C,pars.L3dep,xi,ML));
// 		},
// 		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
//   return phys_res;
// }

// dboot_t cont_chir_fit_m_rat(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_ml_guess={0,0.001};
//   const ave_err_t KPi_guess={0.000,0.001};
//   const ave_err_t K2Pi_guess={-0.000,0.001};
//   const ave_err_t L3dep_guess={0,0.001};
//   const string yaxis_title="$$m_{rat}-1";
  
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
  
//   return cont_chir_fit_m_rat(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,yaxis_title);
// }

// ////////////////////////////////////////////////////////////////////////////////////////////////////////

// template <class Tpars> Tpars FSE_dM(const Tpars &C,const Tpars &L3dep,double L)
// {return C*L3dep/L/L/L/L;}

// //! ansatz fit
// template <class Tpars,class Tm,class Ta>
// Tpars cont_chir_ansatz_dM(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t &isyst)
// {
//   Tpars xi=xi_fun(B0,aml,f0);
  
//   return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_dM(C,L3dep,L);
// }

// //! perform the fit to the continuum limit of mass correction
// dboot_t cont_chir_fit_dM(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &K2Pi_guess,const ave_err_t &L3dep_guess,const string &yaxis_title)
// {
//   //set_printlevel(3);
  
//   boot_fit_t boot_fit;
//   size_t nbeta=a.size();
//   cont_chir_fit_pars_t pars(nbeta);
  
//   pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
//   pars.iC=boot_fit.add_fit_par(pars.C,"C_dM",C_guess.ave(),C_guess.err());
//   pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
//   pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
//   pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  
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
// 		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
//   return phys_res;
// }

// dboot_t cont_chir_fit_dM(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_ml_guess={0,0.001};
//   const ave_err_t KPi_guess={0.000,0.001};
//   const ave_err_t K2Pi_guess={-0.000,0.001};
//   const ave_err_t L3dep_guess={0,0.001};
//   const string yaxis_title="dm";
  
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
  
//   return cont_chir_fit_dM(a,z,f0,B0,ext_data,ml_phys,path,isyst,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,K2Pi_guess,L3dep_guess,yaxis_title);
// }

///////////////////////////////////////////////////////////////////////////////////////////////

class syst_t
{
public:
  double ave,wave,stat,wstat,totsy,tot;
  vector<double> syst;
};

syst_t perform_analysis(const dbvec_t &data,const index_t &ind)
{
  syst_t r;
  double &ave=r.ave=0,&err=r.stat=0,&wave=r.wave=0,&werr=r.wstat=0,weight=0;
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
  
  r.syst=syst_analysis_sep_bis(data.ave_err(),ind);
  double &totsy=r.totsy=0,&tot=r.tot=0;
  for(auto &s : r.syst) totsy+=s*s;
  totsy=sqrt(totsy);
  tot=sqrt(sqr(totsy)+sqr(err));
  
  return r;
}

//! prepare the table for a given quantity
void prepare_table(const string &tag,const dbvec_t &quantity,const index_t &ind,double fact)
{
  index_t ind_table({{"Input",ninput_an},{"HLP",ncont_extrap},{"Fran",nfit_range_variations},{"Inte",nint_num_variations}});
  string path=combine("%s/tables/%s.txt",qname[im].c_str(),tag.c_str());
  ofstream out(path);
  if(!(out.good())) CRASH("Opening %s",path.c_str());
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      dbvec_t data(ind_table.max());
      for(size_t itable=0;itable<ind_table.max();itable++)
	{
	  vector<size_t> comp=ind_table(itable);
	  size_t input_an_id=comp[0],icont_extrap=comp[1],ifit_range=comp[2],iint=comp[3];
	  size_t iintegr=ind({input_an_id,icont_extrap,ifit_range,iens,iint});
	  data[itable]=quantity[iintegr];
	}
      
      syst_t r=perform_analysis(data,ind_table);
      double ave=r.ave*fact;
      vector<double> systs;
      systs.push_back(r.syst[0]*fact);
      systs.push_back(r.syst[1]*fact);
      systs.push_back(sqrt(sqr(r.syst[2])+sqr(r.syst[3]))*fact);
      
      //out<<ens_data[iens].path<<" "<<ave<<endl<<systs<<endl;
      out<<ens_data[iens].path<<" "<<smart_print(ave,systs)<<endl;
    }
}

//! write average, errors and systematics
void perform_analysis(const dbvec_t &data,const index_t &ind,const string &name)
{
  cout<<" ================== "<<name<<" ================== "<<endl;
  auto r=perform_analysis(data,ind);

  for(auto el : vector<pair<string,double>>{{" Ave",r.ave},{"Weig",r.wave}})
    cout<<el.first.substr(0,6)<<":\t"<<el.second<<endl;
  
  vector<pair<string,double>> top;
  top.push_back({" Stat",r.stat});
  top.push_back({"Werr",r.wstat});
  for(size_t i=0;i<ind.rank();i++) top.push_back({ind.name(i),r.syst[i]});
  top.push_back({"TotSy",r.totsy});
  top.push_back({" Tot",r.tot});
  
  for(auto t : top)
    {
      cout<<t.first.substr(0,6);
      cout<<"\t=\t";
      cout<<t.second;
      cout<<"\t=\t";
      streamsize prec=cout.precision();
      cout.precision(3);
      cout<<fabs(t.second/r.ave)*100<<" %";
      cout.precision(prec);
      cout<<endl;
    }
  cout<<endl;
}

int main(int narg,char **arg)
{
  int start=time(0);
  
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
  
  string res_path=combine("%s/integr_results",qname[im].c_str());
  
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
      
      ens_data_t &ens=ens_data[iens];
      size_t TH=ens.T/2;
      string ens_qpath=ens.path+"/plots_"+qname[im];
      
      if(an_mode==compute_everything)
	{
	  djack_t jZ2_V,jZ2_P;
	  two_pts_with_ins_ratio_fit(jZ2_V,jM_V[ind],jA_V[ind],jSL_V[ind],jVV_LO[iens],jVV_QED[iens],TH,tmin_fit(iens,ifit_range),ens.tmax[im],
				     combine("%s/VV_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/VV_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	  two_pts_with_ins_ratio_fit(jZ2_P,jM_P[ind],jA_P[ind],jSL_P[ind],jPP_LO[iens],jPP_QED[iens],TH,tmin_fit(iens,ifit_range),ens.tmax[im],
				     combine("%s/PP_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/PP_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	  jZ_V[ind]=sqrt(jZ2_V);
	  jZ_P[ind]=sqrt(jZ2_P);
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
  index_t ind_integr({{"Input",ninput_an},{"Cont",ncont_extrap},{"FitRange",nfit_range_variations},{"Ens",nens_used},{"IntNum",nint_num_variations}});
  size_t nintegr_max=ind_integr.max();
  dbvec_t resc_a(nintegr_max);
  dbvec_t  LO_correl(nintegr_max), LO_remaind(nintegr_max), LO(nintegr_max);
  dbvec_t QED_correl(nintegr_max),QED_remaind(nintegr_max),QED(nintegr_max);
  dbvec_t                                                  RAT(nintegr_max);
  
  //read if avail
  for(auto &obj : {&resc_a,&LO_correl,&LO_remaind,&QED_correl,&QED_remaind})
    if(an_mode!=compute_everything) obj->bin_read(results_in);
  
  {
    djvec_t A40_QED_data;
    A40_QED_data.resize(3);
    double A40_x[3];
    size_t iA40_L=0;
    for(size_t iens=0;iens<nens_used;iens++)
      {
	ens_data_t &ens=ens_data[iens];
	if(fabs(ens.aml-0.0040)<0.001 and ens.ib==0)
	  {
	    size_t i2pts=ind_2pts_fit({0,iens});
	    double Za=Za_ae[0][0].ave();
	    djack_t a=jM_V[i2pts]/M_V_phys[im];
	    size_t upto=ens.T/2-DT;
	    djack_t c1_QED=integrate_corr_times_kern_up_to(jVV_QED[iens],ens.T,a,im,upto)*sqr(Za);
	    djack_t c2_QED=integrate_QED_reco_from(jA_V[i2pts],jZ_V[i2pts],jSL_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	    A40_QED_data[iA40_L]=c1_QED+c2_QED;
	    
	    A40_x[iA40_L]=1.0/ens.L;
	    
	    iA40_L++;
	  }
      }
    
    {
      string path=combine("%s/tables/four.txt",qname[im].c_str());
      ofstream tab_four(path);
      if(!(tab_four.good())) CRASH("Opening %s",path.c_str());
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  ens_data_t &ens=ens_data[iens];
	  
	  size_t ib=ens.ib;
	  double Za=Za_ae[0][ib].ave();
	  djack_t a;a=1/lat_par[0].ainv[ib].ave();
	  size_t i2pts=ind_2pts_fit({0,iens});
	  size_t tmin=tmin_fit(iens,0),tmax=ens.tmax[im];
	  vector<size_t> possupto({tmin+2,(tmin+tmax)/2,tmax-2,ens.T/2-DT});
	  for(size_t iint=0;iint<nint_num_variations;iint++)
	    {
	      size_t upto=possupto[iint];
	      djack_t c1_LO=integrate_corr_times_kern_up_to(jVV_LO[iens],ens.T,a,im,upto)*sqr(Za);
	      djack_t c2_LO=integrate_LO_reco_from(jZ_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	      djack_t c_LO=c1_LO+c2_LO;
	      tab_four<<ens.path<<"\t"<<upto<<"\t"<<c1_LO.ave_err()<<"\t"<<c2_LO.ave_err()<<"\t"<<c_LO.ave_err()<<endl;
	    }
	  tab_four<<endl;
	  for(size_t iint=0;iint<nint_num_variations;iint++)
	    {
	      size_t upto=possupto[iint];
	      djack_t c1_QED=integrate_corr_times_kern_up_to(jVV_QED[iens],ens.T,a,im,upto)*sqr(Za);
	      djack_t c2_QED=integrate_QED_reco_from(jA_V[i2pts],jZ_V[i2pts],jSL_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	      djack_t c_QED=c1_QED+c2_QED;
	      tab_four<<ens.path<<"\t"<<upto<<"\t"<<c1_QED.ave_err()<<"\t"<<c2_QED.ave_err()<<"\t"<<c_QED.ave_err()<<endl;
	    }
	  tab_four<<"============================================================================"<<endl;
	}
      }
    
    djvec_t C_QED(3),C_dMV(3);
    {
      // djvec_t &C=C_QED;
      double Xa=1/A40_x[0],Xb=1/A40_x[1],Xc=1/A40_x[2];
      djack_t Ya=A40_QED_data[0],Yb=A40_QED_data[1],Yc=A40_QED_data[2];
      
      djack_t Y=(Ya-Yb)/(Ya-Yc);
      grace_file_t fuff(combine("%s/plots/fuff.xmg",qname[im].c_str()));
      fuff.write_line([Xa,Xb,Xc](double x){return (pow(Xa,x)-pow(Xb,x))/(pow(Xa,x)-pow(Xc,x));},-36.05,36.05);
      fuff.new_data_set();
      fuff.write_constant_band(-36.05,36.05,Y);
      
      string path_triplet=combine("%s/triplet.txt",qname[im].c_str());
      ofstream triplet(path_triplet);
      if(!(triplet.good())) CRASH("Opening %s",path_triplet.c_str());
      for(size_t i=0;i<njacks;i++)
	{
	  triplet<<Ya[i]<<endl;
	  triplet<<Yb[i]<<endl;
	  triplet<<Yc[i]<<endl;
	  triplet<<endl;
	}
      // C[2]=Brent_solve<djack_t>([&Y,Xa,Xb,Xc](double nu,size_t ijack){return Y[ijack]-(pow(Xa,nu)-pow(Xb,nu))/(pow(Xa,nu)-pow(Xc,nu));},-36.05,36.05);
      // C[1]=(Ya-Yb)/(pow(Xa,C[2])-pow(Xb,C[2]));
      // C[0]=Ya-C[1]*pow(Xa,C[2]);
    }
    
    cout<<"Exponent QED: "<<C_QED[2].ave_err()<<endl;
    cout<<"Exponent dMV: "<<C_dMV[2].ave_err()<<endl;
    
    grace_file_t fit_file(combine("%s/plots/A40_QED.xmg",qname[im].c_str()));
    fit_file.set_xaxis_label("1/L");
    
    //band of the fit
    fit_file.write_polygon([&C_QED](double x) -> dboot_t {return C_QED[0]+C_QED[1]*pow(x,C_QED[2]);},0.0001,0.055);
    fit_file.new_data_set();
    fit_file.set_settype(grace::XYDY);
    fit_file.new_data_set();
    for(size_t iL=0;iL<3;iL++) fit_file<<A40_x[iL]<<" "<<A40_QED_data[iL].ave_err()<<endl;
  }
  
  const double rat_perturb=-0.01835*sqr(eq[im]);
  for(size_t ind=0;ind<nintegr_max;ind++)
    {
      vector<size_t> comp=ind_integr(ind);
      size_t input_an_id=comp[0],icont_extrap=comp[1],ifit_range=comp[2],iens=comp[3],iint=comp[4];
      ens_data_t &ens=ens_data[iens];
      string ens_qpath=ens.path+"/plots_"+qname[im];
      size_t ib=ens.ib,TH=ens.T/2;
      write_ens_header(iens);
      cout<<"Flags: input_an_id="<<input_an_id<<", icont_extrap="<<icont_extrap<<", ifit_range="<<ifit_range<<", iens="<<iens<<", iint="<<iint<<endl;
      
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
	  
	  //decide the point where to stop the numerical integration
	  size_t tmin=tmin_fit(iens,ifit_range),tmax=ens.tmax[im];
	  size_t upto=vector<size_t>({tmin+2,(tmin+tmax)/2,tmax-2,TH-DT})[iint];
	  
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
      RAT[ind]=QED[ind]/LO[ind]+rat_perturb;
      cout<<" Ratio: "<<RAT[ind].ave_err()<<endl;
	  
      
      cout<<"amu: "<< LO_correl[ind].ave_err()<<" + "<<LO_remaind[ind].ave_err()<<" = "<<LO[ind].ave_err()<<endl;
      cout<<"amu_QED: "<<QED_correl[ind].ave_err()<<" + "<<QED_remaind[ind].ave_err()<<" = "<<QED[ind].ave_err()<<endl;
    }
  
  //write if computed
  for(auto &obj : {&resc_a,&LO_correl,&LO_remaind,&QED_correl,&QED_remaind})
    if(an_mode==compute_everything) obj->bin_write(results_out);
  results_out.close();
  
  prepare_table("LO",LO,ind_integr,1e10);
  prepare_table("QED",QED,ind_integr,1e12);
  prepare_table("RAT",RAT,ind_integr,1e5);
  
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
      // size_t chir_an_id=case_of<c_chir>(isyst);
      // size_t FSE_an_id=case_of<c_FSE>(isyst);
      size_t cont_an_id=case_of<c_cont>(isyst);
      size_t fit_range_an_id=case_of<c_fit_range>(isyst);
      size_t int_num_an_id=case_of<c_int_num>(isyst);
      for(size_t i=0;i<ind_syst.rank();i++) cout<<ind_syst.name(i)<<":\t"<<comp[i]<<endl;
      
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
	{
	  ens_data_t &ens=ens_data[iens];
	  string ens_qpath=ens.path+"/plots_"+qname[im];
	  
	  size_t iintegr=ind_integr({input_an_id,cont_an_id,fit_range_an_id,iens,int_num_an_id});
	  
	  //passing dummy aux mass
	  dboot_t dum;
	  dum=0.0;
	  
	  if(case_of<c_FSE>(isyst)%2!=0 or ens_data[iens].use_for_L)
	    data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,LO[iintegr],LO[iintegr]));
	  if(case_of<c_FSE>(isyst) or ens_data[iens].use_for_L)
	    {
	      data_QED.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,QED[iintegr],QED[iintegr]));
	      data_RAT.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,RAT[iintegr],RAT[iintegr]));
	      
	      // //cc ratio
	      // dboot_t m_rat=M_V/M_P[iens]-1.0;
	      // data_m_rat.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,m_rat,m_rat));
	      
	      // //mass corrections
	      // data_dMP.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,SL_P,SL_P));
	      // data_dMV.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,dum,ens_data[iens].ib,ens_data[iens].L,SL_V,SL_V));
	    }
	}
      
      cLO[isyst]=cont_chir_fit_LO(alist,zlist,f0,B0,data_LO,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_LO_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      cQED[isyst]=cont_chir_fit_QED(alist,zlist,f0,B0,data_QED,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_QED_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      cRAT[isyst]=cont_chir_fit_RAT(alist,zlist,f0,B0,data_RAT,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_RAT_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      
      // c_m_rat[iai]=cont_chir_fit_m_rat(alist,zlist,f0,B0,data_m_rat,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_m_rat_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
      
      // c_dMP[iai]=cont_chir_fit_dM(alist,zlist,f0,B0,data_dMP,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_dMP_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
      // c_dMV[iai]=cont_chir_fit_dM(alist,zlist,f0,B0,data_dMV,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_dMV_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
    }
  
  perform_analysis(cLO,ind_syst,"LO");
  perform_analysis(cQED,ind_syst,"QED");
  perform_analysis(cRAT,ind_syst,"RAT");
  
  // perform_analysis(c_m_rat,"m_rat");
  
  cout<<endl<<"Total time: "<<time(0)-start<<" s to perform "<<ind_syst.max()*(nboots+1)<<" fits"<<endl;
  
  return 0;
}
