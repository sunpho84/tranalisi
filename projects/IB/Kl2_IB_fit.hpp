#ifndef _KL2_IB_FIT_HPP
#define _KL2_IB_FIT_HPP

#include <fit.hpp>
#include <oper.hpp>

const size_t nan_syst=8;
const size_t asyst_mask=1;
const size_t fsyst_mask=2;
const size_t csyst_mask=4;

inline bool FSE_an(size_t an_flag)
{return an_flag & fsyst_mask;}

inline bool cont_an(size_t an_flag)
{return an_flag & asyst_mask;}

inline bool chir_an(size_t an_flag)
{return an_flag & csyst_mask;}

class cont_chir_fit_data_t
{
public:
  double aml,ams;
  size_t ib,L;
  dboot_t wfse,wofse;
  cont_chir_fit_data_t(double aml,double ams,size_t ib,size_t L,dboot_t wfse,dboot_t wofse) : aml(aml),ams(ams),ib(ib),L(L),wfse(wfse),wofse(wofse) {}
};

dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag);

dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag);

dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag);

dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag);

dboot_t cont_chir_quad_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag);

dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag);

dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag);

#endif
