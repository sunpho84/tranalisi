#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr

enum an_mode_t{compute_everything,read_integrals};
an_mode_t an_mode=compute_everything;

bool is_A40(double aml,size_t ib)
{return fabs(aml-0.0040)<0.001 and ib==0;}

template <class T> bool is_A40(const T &ens)
{return is_A40(ens.aml,ens.ib);}

template <class Tpars> Tpars FSE_LO(const Tpars &C,const Tpars &L3dep,const Tpars &xi,const Tpars &M,size_t L,size_t FSE_case)
{
  switch(FSE_case)
    {
    case 0: return C*0.0; break;
    case 1: return C*L3dep*xi*exp(-M*L)/(M*L); break;
    case 2: return C*L3dep/L/L; break;
    default: CRASH("Unknown case");return C;break;
    }
}

//! return the tmin for the fit
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

//! approximated Za_fact
dboot_t Za_fact;
extern djack_t *Zm_fact;

//! return the perturbative correction to Za
dboot_t Za_perturb_QED(size_t im)
{return -0.01835*sqr(eq[im])*include_ZA_perturb*Za_fact;}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_LO(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t &isyst)
{
  Tpars xi=xi_fun(B0,aml,f0),M=M_fun(B0,aml);
  return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+xi*adep_ml))+FSE_LO(C,L3dep,xi,M,L,case_of<c_FSE>(isyst));
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
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
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
		[&ext_data,&pars,&B0,&f0,isyst]
		(size_t idata,bool without_with_fse,size_t ib)
		{
		  dboot_t aml=ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib];
		  dboot_t xi=xi_fun(pars.fit_B0,aml,pars.fit_f0),M=M_fun(pars.fit_B0,aml);
		  size_t L=ext_data[idata].L;
		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_LO(pars.C,pars.L3dep,xi,M,L,case_of<c_FSE>(isyst)));},
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
  switch(case_of<c_FSE>(isyst))
    {
    case 0:boot_fit.fix_par_to(pars.iL4dep,0);break;
    case 1:boot_fit.fix_par_to(pars.iL4dep,-2);break;
    case 2:boot_fit.fix_par_to(pars.iL4dep,-1);break;
    }
  
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],p[pars.iL4dep]);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,pars.L4dep);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(combine(path.c_str(),"cont_chir"),ext_data,pars,
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
  
  grace_file_t out_FSE(combine(path.c_str(),"FSE"));
  out_FSE.set_subtitle(ind_syst.descr(isyst));
  out_FSE.set_color_scheme({grace::RED,grace::BLUE});
  out_FSE.set_symbol_scheme({grace::SQUARE,grace::CIRCLE});
  double mlA40=dboot_t(0.0040/pars.fit_a[0]/pars.fit_z[0]).ave(),aA=pars.fit_a[0].ave();
  for(size_t without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      out_FSE.write_polygon([&pars,mlA40,aA,without_with_fse](double x)
    {return cont_chir_ansatz_QED(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,mlA40,aA,pars.adep,pars.adep_ml,1.0/x,pars.L3dep*
				 without_with_fse,pars.L4dep);},1e-5,1.0/12);
      out_FSE.new_data_set();
      
      out_FSE.set_settype(grace::XYDY);
      for(size_t iens=0;iens<nens_used;iens++)
	 if(is_A40(ext_data[iens]))
	   out_FSE<<1.0/ext_data[iens].L<<" "<<dboot_t(ext_data[iens].wfse-without_with_fse*FSE_QED(pars.C,pars.L3dep,ext_data[iens].L,pars.L4dep)).ave_err()<<endl;
      out_FSE.new_data_set();
    }
  
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
      C_guess=ave_err_t({-3e-12,0.5e-12});
      adep_guess={-1.5,0.6};
      KPi_guess={0.0,0.01};
      break;
    case icharm:
      C_guess=ave_err_t({0-4e-12,2e-12});
      adep_guess={-3,0.2};
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
      C_guess=ave_err_t({-0.0005,0.0004});
      adep_guess={-1,0.2};
  break;
    case icharm:
      C_guess=ave_err_t({-0.003,0.002});
      adep_guess={-3,0.3};
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
      systs.push_back(r.stat*fact);
      systs.push_back(r.syst[0]*fact);
      systs.push_back(r.syst[1]*fact);
      systs.push_back(sqrt(sqr(r.syst[2])+sqr(r.syst[3]))*fact);
      
      //out<<ens_data[iens].path<<" "<<ave<<endl<<systs<<endl;
      out<<ens_data[iens].path<<" "<<smart_print(ave,systs)<<endl;
    }
}

//! prepare the table for a given quantity, with or without hl
void prepare_table_with_without_hl(const string &tag,const dbvec_t &quantity,const index_t &ind,double fact,size_t icont_extrap)
{
  index_t ind_table({{"Input",ninput_an},{"Fran",nfit_range_variations},{"Inte",nint_num_variations}});
  string path=combine("%s/tables/%s_%s_hl.txt",qname[im].c_str(),tag.c_str(),(icont_extrap==0)?"without":"with");
  ofstream out(path);
  if(!(out.good())) CRASH("Opening %s",path.c_str());
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      dbvec_t data(ind_table.max());
      for(size_t itable=0;itable<ind_table.max();itable++)
	{
	  vector<size_t> comp=ind_table(itable);
	  size_t input_an_id=comp[0],ifit_range=comp[1],iint=comp[2];
	  size_t iintegr=ind({input_an_id,icont_extrap,ifit_range,iens,iint});
	  data[itable]=quantity[iintegr];
	}
      
      syst_t r=perform_analysis(data,ind_table);
      double ave=r.ave*fact;
      vector<double> systs;
      systs.push_back(r.stat*fact);
      systs.push_back(r.syst[0]*fact);
      systs.push_back(r.syst[1]*fact);
      systs.push_back(sqrt(sqr(r.syst[2])+sqr(r.syst[3]))*fact);
      
      //out<<ens_data[iens].path<<" "<<smart_print(ave,systs)<<endl;
      out<<ens_data[iens].path<<" "<<ave<<" "<<systs[0]<<" "<<systs[1]<<" "<<systs[2]<<endl;
    }
}


//! test the exponents of the A40
void test_exponent_A40(const index_t &ind_2pts_fit,const djvec_t &jM_P,const djvec_t &jM_V,const vector<djvec_t> &jVV_LO,const vector<djvec_t> &jVV_QED,
		       const djvec_t &k_DZ2_rel_V,const djvec_t &jZ2_V,const djvec_t &jSL_V)
{
  size_t nA40=count_if(ens_data.begin(),ens_data.end(),is_A40<ens_data_t>);
  size_t nA=count_if(ens_data.begin(),ens_data.end(),[](const ens_data_t &ens){return ens.ib==0;});
  
  djvec_t A_QED_data(nA);
  djvec_t A40_MP_data(nA40);
  djvec_t A40_MV_data(nA40);
  djvec_t A40_QED_data(nA40);
  djvec_t A40_LO_data(nA40);
  double A_x[nA];
  double A40_x[nA40];
  size_t iA_L=0,iA40_L=0;
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      size_t ib=ens.ib;
      if(ib==0)
	{
	  size_t i2pts=ind_2pts_fit({0,iens});
	  double Za=Za_ae[0][ib].ave();
	  djack_t a=djack_t(jM_V[i2pts]/M_V_phys[im]);
	  cout<<"a for ensemble "<<ens.path<<": "<<smart_print(a.ave_err())<<endl;
	  size_t upto=tmin_fit(iens,0);
	  djack_t c1_QED=integrate_corr_times_kern_up_to(jVV_QED[iens],ens.T,a,im,upto)*sqr(Za);
	  djack_t c2_QED=integrate_QED_reco_from(k_DZ2_rel_V[i2pts],jZ2_V[i2pts],jSL_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	  djack_t c1_LO=integrate_corr_times_kern_up_to(jVV_LO[iens],ens.T,a,im,upto)*sqr(Za);
	  djack_t c2_LO=integrate_LO_reco_from(jZ2_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	  
	  {
	    A_QED_data[iA_L]=(c1_QED+c2_QED)+Za_perturb_QED(im)*(c1_LO+c2_LO);
	    A_x[iA_L]=1.0/ens.L+ens.aml/10;
	    iA_L++;
	  }
	  
	  if(is_A40(ens))
	    {
	      A40_MP_data[iA40_L]=jM_P[i2pts];
	      A40_MV_data[iA40_L]=jM_V[i2pts];
	      
	      A40_LO_data[iA40_L]=c1_LO+c2_LO;
	      A40_QED_data[iA40_L]=c1_QED+c2_QED+Za_perturb_QED(im)*A40_LO_data[iA40_L];
	      
	      A40_x[iA40_L]=1.0/ens.L;
	      
	      iA40_L++;
	    }
	}
    }
  
  djvec_t C_QED_A40(3);
  {
    jack_fit_t fitter;
    C_QED_A40[1]=(A40_QED_data[nA40-1]-A40_QED_data[1])/(A40_x[nA40-1]-A40_x[1]);
    C_QED_A40[0]=A40_QED_data[1]-A40_x[1]*C_QED_A40[1];
    C_QED_A40[2].fill_gauss(1,1.0,624134);
    
    cout<<"Guess: "<<endl;
    cout<<C_QED_A40.ave_err()<<endl;
    
    fitter.add_fit_par(C_QED_A40[0],"C0",C_QED_A40[0].ave_err());
    fitter.add_fit_par(C_QED_A40[1],"C1",C_QED_A40[1].ave_err());
    fitter.add_fit_par(C_QED_A40[2],"C2",{1.0,1.0});
    for(size_t i=0;i<nA40;i++)
      fitter.add_point(A40_QED_data[i],
		       [&A40_x,i](const vector<double> &p,int iel)
		       {return p[0]+p[1]*pow(A40_x[i],p[2]);});
    
    fitter.fit();
  }
  cout<<"Exponent QED from A40: "<<C_QED_A40[2].ave_err()<<endl;
  
  djvec_t C_QED_A(3);
  {
    jack_fit_t fitter;
    C_QED_A[1]=(A_QED_data[nA-1]-A_QED_data[0])/(A_x[nA-1]-A_x[0]);
    C_QED_A[0]=A_QED_data[0]-A_x[0]*C_QED_A[1];
    C_QED_A[2].fill_gauss(1,1.0,624134);
    
    cout<<"Guess: "<<endl;
    cout<<C_QED_A.ave_err()<<endl;
    
    fitter.add_fit_par(C_QED_A[0],"C0",C_QED_A[0].ave_err());
    fitter.add_fit_par(C_QED_A[1],"C1",C_QED_A[1].ave_err());
    fitter.add_fit_par(C_QED_A[2],"C2",C_QED_A[2].ave_err());
    for(size_t i=0;i<nA;i++)
      fitter.add_point(A_QED_data[i],
		       [&A_x,i](const vector<double> &p,int iel)
		       {return p[0]+p[1]*pow(A_x[i],p[2]);});
    
    fitter.fit();
  }
  
  cout<<C_QED_A.ave_err()<<endl;
  
  cout<<"Exponent QED from all A: "<<C_QED_A[2].ave_err()<<endl;
  
  grace_file_t fit_file_QED_A40(combine("%s/plots/A40_QED.xmg",qname[im].c_str()));
  fit_file_QED_A40.set_xaxis_label("1/L");
  fit_file_QED_A40.write_polygon([&C_QED_A40](double x) -> dboot_t {return C_QED_A40[0]+C_QED_A40[1]*pow(x,C_QED_A40[2]);},0.0001,*max_element(A40_x,A40_x+nA40)*1.1);
  fit_file_QED_A40.new_data_set();
  fit_file_QED_A40.set_settype(grace::XYDY);
  fit_file_QED_A40.new_data_set();
  for(size_t iL=0;iL<nA40;iL++) fit_file_QED_A40<<A40_x[iL]<<" "<<A40_QED_data[iL].ave_err()<<endl;
  
  grace_file_t fit_file_QED_A(combine("%s/plots/A_QED.xmg",qname[im].c_str()));
  fit_file_QED_A.set_xaxis_label("1/L");
  fit_file_QED_A.write_polygon([&C_QED_A](double x) -> dboot_t {return C_QED_A[0]+C_QED_A[1]*pow(x,C_QED_A[2]);},0.0001,*max_element(A_x,A_x+nA)*1.1);
  fit_file_QED_A.new_data_set();
  fit_file_QED_A.set_settype(grace::XYDY);
  fit_file_QED_A.new_data_set();
  for(size_t iL=0;iL<nA;iL++) fit_file_QED_A<<A_x[iL]<<" "<<A_QED_data[iL].ave_err()<<endl;
  
  grace_file_t fit_file_LO(combine("%s/plots/A40_LO.xmg",qname[im].c_str()));
  fit_file_LO.set_xaxis_label("1/L");
  fit_file_LO.set_settype(grace::XYDY);
  for(size_t iL=0;iL<nA40;iL++) fit_file_LO<<A40_x[iL]<<" "<<A40_LO_data[iL].ave_err()<<endl;
  
  grace_file_t fit_file_RAT(combine("%s/plots/A40_RAT.xmg",qname[im].c_str()));
  fit_file_RAT.set_xaxis_label("1/L");
  fit_file_RAT.set_settype(grace::XYDY);
  for(size_t iL=0;iL<nA40;iL++) fit_file_RAT<<A40_x[iL]<<" "<<djack_t(A40_QED_data[iL]/A40_LO_data[iL]).ave_err()<<endl;
  
  grace_file_t fit_file_MP(combine("%s/plots/A40_MP.xmg",qname[im].c_str()));
  fit_file_MP.set_xaxis_label("1/L");
  fit_file_MP.set_settype(grace::XYDY);
  for(size_t iL=0;iL<nA40;iL++) fit_file_MP<<A40_x[iL]<<" "<<A40_MP_data[iL].ave_err()<<endl;
  
  grace_file_t fit_file_MV(combine("%s/plots/A40_MV.xmg",qname[im].c_str()));
  fit_file_MV.set_xaxis_label("1/L");
  fit_file_MV.set_settype(grace::XYDY);
  for(size_t iL=0;iL<nA40;iL++) fit_file_MV<<A40_x[iL]<<" "<<A40_MV_data[iL].ave_err()<<endl;
}

int main(int narg,char **arg)
{
  int start=time(0);
  
  gm2_initialize(narg,arg);
  Za_fact.fill_gauss({0.9,0.1,9873834});
  Zm_fact=new djack_t;
  Zm_fact->fill_gauss({1.0,1e-6,6232342});
  
  vector<djvec_t> jPP_LO(nens_used),jPP_QED(nens_used),jVV_LO(nens_used),jVV_QED(nens_used);
  string deltam_cr_path=qname[im]+"/tables/deltam_cr.txt";
  ofstream deltam_cr_file(deltam_cr_path);
  if(not deltam_cr_file.good()) CRASH("Unable to open %s",deltam_cr_path.c_str());
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      djack_t deltam_cr;
      deltam_cr=compute_deltam_cr(ens,icharm);
      deltam_cr=compute_deltam_cr(ens,istrange);
      deltam_cr=compute_deltam_cr(ens,ilight);
      deltam_cr.bin_write(ens.path+"/deltam_cr_light.bin");
      deltam_cr_file<<ens.path<<"\t"<<deltam_cr.ave_err()<<endl;
      
      //load LO
      jPP_LO[iens]=read_PP("00",ens,im,1,RE);
      jVV_LO[iens]=read_VV("00",ens,im,1,RE);
      jPP_LO[iens].ave_err().write(ens.path+"/plots_"+qname[im]+"/jPP_LO_corr.xmg");
      jVV_LO[iens].ave_err().write(ens.path+"/plots_"+qname[im]+"/jVV_LO_corr.xmg");
      
      //test numerical
      string numpath=ens.path+"/data/corrNmNm_V1V1_"+qname[im][0]+qname[im][0];
      cout<<"Checking file: "<<numpath<<endl;
      if(file_exists(numpath))
	{
	  double coef=16;
	  djvec_t NmNm=djvec_t(read_VV("NmNm",ens,im,1,RE)-jVV_LO[iens])*coef;
	  djvec_t NpNp=djvec_t(read_VV("NpNp",ens,im,1,RE)-jVV_LO[iens])*coef;
	  NmNm.ave_err().write(ens.path+"/plots_"+qname[im]+"/jVV_QED_num_m.xmg");
	  NpNp.ave_err().write(ens.path+"/plots_"+qname[im]+"/jVV_QED_num_p.xmg");
	  djvec_t(0.5*NmNm+0.5*NpNp).ave_err().write(ens.path+"/plots_"+qname[im]+"/jVV_QED_num.xmg");
	  
	  NmNm=djvec_t(read_PP("NmNm",ens,im,1,RE)-jPP_LO[iens])*coef;
	  NpNp=djvec_t(read_PP("NpNp",ens,im,1,RE)-jPP_LO[iens])*coef;
	  NmNm.ave_err().write(ens.path+"/plots_"+qname[im]+"/jPP_QED_num_m.xmg");
	  NpNp.ave_err().write(ens.path+"/plots_"+qname[im]+"/jPP_QED_num_p.xmg");
	  djvec_t(0.5*NmNm+0.5*NpNp).ave_err().write(ens.path+"/plots_"+qname[im]+"/jPP_QED_num.xmg");
	}
      
      double a=1/lat_par[0].ainv[ens.ib].ave();
      //double a=constant_fit(effective_mass(jVV_LO[iens])/M_V_phys[im],ens.tmin[im],ens.tmax[im],ens.path+"/plots_"+qname[im]+"/test_a_read_QED.xmg").ave();
      cout<<"Using a="<<a<<" to compute QED contribution for ensemble "<<ens.path<<endl;
      
      //load QED
      jPP_QED[iens]=read_QED_PP(ens,1,im,deltam_cr,jPP_LO[iens],a);
      jVV_QED[iens]=read_QED_VV(ens,1,im,deltam_cr,jVV_LO[iens],a);
      jPP_QED[iens].ave_err().write(ens.path+"/plots_"+qname[im]+"/jPP_QED_corr.xmg");
      jVV_QED[iens].ave_err().write(ens.path+"/plots_"+qname[im]+"/jVV_QED_corr.xmg");
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
  djvec_t jZ2_V(n2pts_fit_max),jM_V(n2pts_fit_max),k_DZ2_rel_V(n2pts_fit_max),jSL_V(n2pts_fit_max);
  djvec_t jZ2_P(n2pts_fit_max),jM_P(n2pts_fit_max),k_DZ2_rel_P(n2pts_fit_max),jSL_P(n2pts_fit_max);
  
  //read if avail
  for(auto &obj : {&jZ2_P,&jM_P,&k_DZ2_rel_P,&jSL_P,&jZ2_V,&jM_V,&k_DZ2_rel_V,&jSL_V})
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
	  grace_file_t out(combine("%s/VV_QED_an%zu_jacks.xmg",ens_qpath.c_str(),ifit_range));
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    {
	      for(size_t t=0;t<jVV_QED[iens].size();t++) out<<t<<" "<<jVV_QED[iens][t][njacks]*njacks-jVV_QED[iens][t][ijack]*(njacks-1)<<endl;
	      out<<endl;
	    }
	  
	  two_pts_with_ins_ratio_fit(jZ2_V[ind],jM_V[ind],k_DZ2_rel_V[ind],jSL_V[ind],jVV_LO[iens],jVV_QED[iens],TH,tmin_fit(iens,ifit_range),ens.tmax[im],
				     combine("%s/VV_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/VV_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	  two_pts_with_ins_ratio_fit(jZ2_P[ind],jM_P[ind],k_DZ2_rel_P[ind],jSL_P[ind],jPP_LO[iens],jPP_QED[iens],TH,tmin_fit(iens,ifit_range),ens.tmax[im],
				     combine("%s/PP_LO_an%zu.xmg",ens_qpath.c_str(),ifit_range),combine("%s/PP_QED_an%zu.xmg",ens_qpath.c_str(),ifit_range));
	}
      
      cout<<"aM_V: "<<jM_V[ind].ave_err()<<", adeltaM_V/e2: "<<djack_t(-jSL_V[ind]/e2).ave_err()<<endl;
      cout<<"aM_P: "<<jM_P[ind].ave_err()<<", adeltaM_P/e2: "<<djack_t(-jSL_P[ind]/e2).ave_err()<<endl;
    }
  
  //write if computed
  for(auto &obj : {&jZ2_P,&jM_P,&k_DZ2_rel_P,&jSL_P,&jZ2_V,&jM_V,&k_DZ2_rel_V,&jSL_V})
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
	vector<djack_t> c_LO(nint_num_variations);
	vector<djack_t> c1_LO(nint_num_variations);
	vector<djack_t> c2_LO(nint_num_variations);
	for(size_t iint=0;iint<nint_num_variations;iint++)
	  {
	    size_t upto=possupto[iint];
	    c1_LO[iint]=integrate_corr_times_kern_up_to(jVV_LO[iens],ens.T,a,im,upto)*sqr(Za);
	    c2_LO[iint]=integrate_LO_reco_from(jZ2_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za);
	    c_LO[iint]=c1_LO[iint]+c2_LO[iint];
	    tab_four<<ens.path<<"\t"<<upto<<"\t"<<c1_LO[iint].ave_err()<<"\t"<<c2_LO[iint].ave_err()<<"\t"<<c_LO[iint].ave_err()<<endl;
	  }
	tab_four<<endl;
	for(size_t iint=0;iint<nint_num_variations;iint++)
	  {
	    size_t upto=possupto[iint];
	    djack_t c1_QED=integrate_corr_times_kern_up_to(jVV_QED[iens],ens.T,a,im,upto)*sqr(Za)+c1_LO[iint]*Za_perturb_QED(im).ave();
	    djack_t c2_QED=integrate_QED_reco_from(k_DZ2_rel_V[i2pts],jZ2_V[i2pts],jSL_V[i2pts],jM_V[i2pts],a,im,upto)*sqr(Za)+c2_LO[iint]*Za_perturb_QED(im).ave();
	    djack_t c_QED=c1_QED+c2_QED;
	    tab_four<<ens.path<<"\t"<<upto<<"\t"<<c1_QED.ave_err()<<"\t"<<c2_QED.ave_err()<<"\t"<<c_QED.ave_err()<<endl;
	  }
	tab_four<<"============================================================================"<<endl;
      }
    test_exponent_A40(ind_2pts_fit,jM_P,jM_V,jVV_LO,jVV_QED,k_DZ2_rel_V,jZ2_V,jSL_V);
  }
  
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
      dboot_t Z2_V(bi,jZ2_V[i2pts]),M_V(bi,jM_V[i2pts]),A_V(bi,k_DZ2_rel_V[i2pts]),SL_V(bi,jSL_V[i2pts]);
      
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
	  LO_remaind[ind]=integrate_LO_reco_from(Z2_V,M_V,resc_a[ind],im,upto)*sqr(Za[ib]);
	  
	  QED_correl[ind]=integrate_corr_times_kern_up_to(VV_QED,ens.T,resc_a[ind],im,upto)*sqr(Za[ib])+Za_perturb_QED(im)*LO_correl[ind];
	  QED_remaind[ind]=integrate_QED_reco_from(A_V,Z2_V,SL_V,M_V,resc_a[ind],im,upto)*sqr(Za[ib])+Za_perturb_QED(im)*LO_remaind[ind];
	  
	  index_t ind_rest({{"Input",ninput_an},{"Cont",ncont_extrap},{"FitRange",nfit_range_variations}});
	  size_t irest=ind_rest({input_an_id,icont_extrap,ifit_range});
	  
	  compare_LO_num_reco(combine("%s/kern_LO_num_reco_%zu.xmg",ens_qpath.c_str(),irest),VV_LO,Z2_V,M_V,resc_a[iens]);
	  compare_QED_num_reco(combine("%s/kern_QED_num_reco_%zu.xmg",ens_qpath.c_str(),irest),VV_QED,A_V,Z2_V,SL_V,M_V,resc_a[iens]);
	}
      
      LO[ind]=LO_correl[ind]+LO_remaind[ind];
      QED[ind]=QED_correl[ind]+QED_remaind[ind];
      RAT[ind]=QED[ind]/LO[ind];
      cout<<" Ratio: "<<RAT[ind].ave_err()<<endl;
      
      cout<<"amu: "<< LO_correl[ind].ave_err()<<" + "<<LO_remaind[ind].ave_err()<<" = "<<LO[ind].ave_err()<<endl;
      cout<<"amu_QED (with QED perturb correction to Za): "<<QED_correl[ind].ave_err()<<" + "<<QED_remaind[ind].ave_err()<<" = "<<QED[ind].ave_err()<<endl;
    }
  
  //write if computed
  for(auto &obj : {&resc_a,&LO_correl,&LO_remaind,&QED_correl,&QED_remaind})
    if(an_mode==compute_everything) obj->bin_write(results_out);
  results_out.close();

  const double qed_print_fact[3]={1e12,1e12,1e13};
  prepare_table("LO",LO,ind_integr,1e10);
  prepare_table("QED",QED,ind_integr,qed_print_fact[im]);
  prepare_table("RAT",RAT,ind_integr,1e5);
  for(size_t icont_extrap=0;icont_extrap<2;icont_extrap++)
    {
      prepare_table_with_without_hl("LO",LO,ind_integr,1e10,icont_extrap);
      prepare_table_with_without_hl("QED",QED,ind_integr,qed_print_fact[im],icont_extrap);
      prepare_table_with_without_hl("RAT",RAT,ind_integr,1e5,icont_extrap);
    }
  
  
  //loop over analysis flags and input scale determination
  dbvec_t cLO(ind_syst.max());
  dbvec_t cQED(ind_syst.max());
  dbvec_t cRAT(ind_syst.max());
  //dbvec_t c_m_rat,c_dMP;
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
	  grace::default_color_scheme={grace::RED,grace::RED,grace::RED,grace::RED, grace::BLUE,grace::BLUE, grace::GREEN4,grace::VIOLET};
	  grace::default_symbol_scheme={grace::STAR,grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
	  // grace::default_color_scheme={grace::RED,grace::RED,grace::RED, grace::BLUE,grace::BLUE, grace::GREEN4,grace::VIOLET};
	  // grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
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
	  
	  if(case_of<c_FSE>(isyst)!=0 or ens_data[iens].use_for_L)
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
      cQED[isyst]=cont_chir_fit_QED(alist,zlist,f0,B0,data_QED,lat_par[input_an_id].ml,combine("%s/plots/%s_QED_an%zu.xmg",qname[im].c_str(),"%s",isyst),isyst,false,beta_list);
      cRAT[isyst]=cont_chir_fit_RAT(alist,zlist,f0,B0,data_RAT,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_RAT_an%zu.xmg",qname[im].c_str(),isyst),isyst,false,beta_list);
      
      // c_m_rat[iai]=cont_chir_fit_m_rat(alist,zlist,f0,B0,data_m_rat,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_m_rat_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
      
      // c_dMP[iai]=cont_chir_fit_dM(alist,zlist,f0,B0,data_dMP,lat_par[input_an_id].ml,combine("%s/plots/cont_chir_dMP_flag%zu_an%zu.xmg",qname[im].c_str(),isyst,input_an_id),isyst,false,beta_list);
    }
  
  perform_analysis(cLO,ind_syst,"LO");
  perform_analysis(cQED,ind_syst,"QED");
  perform_analysis(cRAT,ind_syst,"RAT");
  
  // perform_analysis(c_m_rat,"m_rat");
  
  cout<<endl<<"Total time: "<<time(0)-start<<" s to perform "<<ind_syst.max()*(nboots+1)<<" fits"<<endl;
  
  return 0;
}
