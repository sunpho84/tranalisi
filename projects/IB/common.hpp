#ifndef _COMMON_HPP
#define _COMMON_HPP

#include <cstdlib>
#include <tranalisi.hpp>

#ifndef EXTERN_COMMON
 #define EXTERN_COMMON extern
 #define INIT_TO(A)
 #define INIT_ARGS(A)
#else
 #define INIT_TO(A) =A
 #define INIT_ARGS(...) (__VA_ARGS__)
#endif

//total number of possible ensemble
const size_t nens_total=15;

//number of input analysis
const size_t ninput_an=8;

//! number of beta
const size_t nbeta=3;

//betas
const vector<string> beta_list={"1.90","1.95","2.10"};

//! number of boostrap
const size_t nboots=100;

//! charge of up quark
const double eu=2.0/3;
//! charge of down quark
const double ed=-1.0/3;
//! charge of strange quark
const double es=ed;
//! charge of charm
const double ec=eu;

//! continuum limit
const double a_cont=1e-5;

//! renormalization constants
const vector<vector<ave_err_t>> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004}},{{0.703,0.002},{0.714,0.002},{0.752,0.002}}});
const vector<vector<ave_err_t>> Zt_ae({{{0.711,0.005},{0.724,0.004},{0.774,0.004}},{{0.700,0.003},{0.711,0.002},{0.767,0.002}}});
const int Za_seed[nbeta]={13124,862464,76753};
const int Zt_seed[nbeta]={5634,917453,324338};

EXTERN_COMMON dbvec_t Za,Zt;

//! boot initializer
EXTERN_COMMON boot_init_t bi;

//! list of a
EXTERN_COMMON dbvec_t alist INIT_ARGS(nbeta),zlist INIT_ARGS(nbeta);

//! infinite volum limit
const double inf_vol=1e10;

//! hold the jacknife index for the given bootstrap
EXTERN_COMMON boot_init_t jack_index[ninput_an][nens_total];

//! symbols for plots
const vector<size_t> symbol={grace::SQUARE,grace::CIRCLE,grace::DIAMOND};
//! colors for plots
const vector<size_t> color={grace::GREEN4,grace::RED,grace::BLUE};

//! class to hold results from mass determination
class lat_par_t
{
public:
  dboot_t ml,ms,mc,r0,f0,B0;
  dbvec_t ainv,Z;
 
  lat_par_t() : ainv(nbeta),Z(nbeta) {}
};

//! results from mass determinations
EXTERN_COMMON vector<lat_par_t> lat_par INIT_ARGS(ninput_an);

//! initialize the common part of the IB analsys
void init_common_IB(string ens_pars);

//! read an ascii bootstrap
dboot_t read_boot(const raw_file_t &file);

//! perform the analysis according to eq.28
ave_err_t eq_28_analysis(const dbvec_t &v);

//! hold data for continuum chiral infvol extrapolation
class cont_chir_fit_data_t
{
public:
  double aml,ams;
  dboot_t aMaux;
  size_t ib,L;
  dboot_t wfse,wofse;
  cont_chir_fit_data_t(double aml,double ams,dboot_t aMaux,size_t ib,size_t L,dboot_t wfse,dboot_t wofse) : aml(aml),ams(ams),aMaux(aMaux),ib(ib),L(L),wfse(wfse),wofse(wofse) {}
};

//! holds index and out pars
class cont_chir_fit_pars_t
{
public:
  dbvec_t ori_a,ori_z;
  dboot_t ori_f0,ori_B0;
  dboot_t ori_ml_phys;
  size_t nbeta;
  vector<size_t> ipara,iparz;
  size_t if0,iB0;
  size_t iadep,iadep_ml,iL2dep,iL3dep,iL4dep,iML4dep;
  size_t iC,iKPi,iKK,iK2Pi,iK2K;
  size_t iml_phys;
  bool fitting_a,fitting_z;
  dbvec_t fit_a,fit_z;
  dboot_t fit_C,fit_f0,fit_B0,fit_L2dep,fit_L3dep,fit_L4dep,fit_ML4dep;
  dboot_t adep,adep_ml,L2dep,L3dep,L4dep,ML4dep;
  dboot_t C,KPi,KK,K2Pi,K2K;
  dboot_t fit_ml_phys;
  
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
    iadep=boot_fit.add_fit_par(adep,"adep",adep_guess.ave(),adep_guess.err());
    iadep_ml=boot_fit.add_fit_par(adep_ml,"adep_ml",adep_ml_guess.ave(),adep_ml_guess.err());
  }
  
  //! add fsedep
  void add_fsedep_pars(const ave_err_t &L3dep_guess,const ave_err_t &L4dep_guess,const ave_err_t &ML4dep_guess,boot_fit_t &boot_fit,const bool flag_prior=true)
  {
    if(flag_prior)
      {
	iL3dep=boot_fit.add_fit_par(L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
	iL4dep=boot_fit.add_fit_par(L4dep,"L4dep",L4dep_guess.ave(),L4dep_guess.err());
	iML4dep=boot_fit.add_fit_par(ML4dep,"ML4dep",ML4dep_guess.ave(),ML4dep_guess.err());
      }
    else
      {
	iL3dep=boot_fit.add_self_fitted_point(fit_L3dep,"L3dep_prior",{3.0,2.0});
	iL4dep=boot_fit.add_self_fitted_point(fit_L4dep,"L4dep_prior",{0.0,0.1});
	iML4dep=boot_fit.add_self_fitted_point(fit_ML4dep,"ML4dep_prior",{0.0,0.1});
      }
  }
  
  //! add all common pars
  void add_common_pars(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
		       const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,boot_fit_t &boot_fit)
  {
    add_az_pars(a,z,boot_fit);
    add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
    //f0 and B0
    ori_f0=f0;
    ori_B0=B0;
    if0=boot_fit.add_self_fitted_point(fit_f0,"f0",f0,DO_NOT_CORRELATE);
    iB0=boot_fit.add_self_fitted_point(fit_B0,"B0",B0,DO_NOT_CORRELATE);
  }
  
  //! add the low energy constants
  void add_LEC_pars(const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &KK_guess,const ave_err_t &K2Pi_guess,const ave_err_t &K2K_guess,boot_fit_t &boot_fit,const bool flag_prior=true)
  {
    if(flag_prior)
      {
	iC=boot_fit.add_fit_par(C,"C",C_guess.ave(),C_guess.err());
	
      }
    else
      {
	iC=boot_fit.add_self_fitted_point(fit_C,"C_prior",{3.54e-5,1.0e-5});
      }
    
    iKPi=boot_fit.add_fit_par(KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
    iKK=boot_fit.add_fit_par(KK,"KK",KK_guess.ave(),KK_guess.err());
    iK2Pi=boot_fit.add_fit_par(K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
    iK2K=boot_fit.add_fit_par(K2K,"K2K",K2K_guess.ave(),K2K_guess.err());
    
  }
  
  //! add the dependency from ml_phys
  void add_ml_phys_par(const dboot_t &ml_phys,boot_fit_t &boot_fit)
  {
    ori_ml_phys=ml_phys;
    iml_phys=boot_fit.add_self_fitted_point(fit_ml_phys,"ml_phys",ml_phys,DO_NOT_CORRELATE);
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
    cout<<"ml_phys: "<<fit_ml_phys.ave_err()<<", orig: "<<ori_ml_phys.ave_err()<<", ratio: "<<dboot_t(fit_ml_phys/ori_ml_phys-1).ave_err()<<endl;
    //
    cout<<"Adep: "<<adep.ave_err()<<endl;
    cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
    cout<<"L2dep: "<<L2dep.ave_err()<<endl;
    cout<<"L3dep: "<<L3dep.ave_err()<<endl;
    cout<<"L4dep: "<<L4dep.ave_err()<<endl;
    cout<<"ML4dep: "<<ML4dep.ave_err()<<endl;
    cout<<"L3dep_prior: "<<fit_L3dep.ave_err()<<endl;
    cout<<"L4dep_prior: "<<fit_L4dep.ave_err()<<endl;
    cout<<"ML4dep_prior: "<<fit_ML4dep.ave_err()<<endl;
  }
  
  //! print the value
  void print_LEC_pars() const
  {
    cout<<"C: "<<C.ave_err()<<endl;
    cout<<"C_prior: "<<fit_C.ave_err()<<endl;
    cout<<"KPi: "<<KPi.ave_err()<<endl;
    cout<<"KK: "<<KK.ave_err()<<endl;
    cout<<"K2Pi: "<<K2Pi.ave_err()<<endl;
    cout<<"K2K: "<<K2K.ave_err()<<endl;
  }
};

//! perform the fit to the continuum limit
void cont_chir_fit_minimize
(const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,boot_fit_t &boot_fit,double apow,double zpow,
 const function<double(const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double Maux,double ac,double L)> &cont_chir_ansatz,bool cov_flag);

//! plot the continuum-chiral extrapolation
void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,const string &subtitle="");

void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,size_t univ_full_sub);

void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,size_t univ_full_sub,size_t FSE_flag);

vector<double> syst_analysis_sep_bis(const vector<ave_err_t> &v,const index_t &fact);

#undef INIT_TO
#undef EXTERN_COMMON

#endif
