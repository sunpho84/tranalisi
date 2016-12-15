#ifndef _KL2_IB_FIT_HPP
#define _KL2_IB_FIT_HPP

#include <fit.hpp>

class cont_chir_fit_data_t_pi
{
public:
  double aml;
  size_t ib,L;
  dboot_t y,fse;;
  cont_chir_fit_data_t_pi(double aml,size_t ib,size_t L,dboot_t y,dboot_t fse) : aml(aml),ib(ib),L(L),y(y),fse(fse) {}
};

//! perform a fit to the continuum and chiral
void cont_chir_fit_pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_pi> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an);

// class cont_chir_fit_data_t_epsilon
// {
// public:
//   double aml,ams;
//   size_t ib,L;
//   dboot_t wfse,wofse;
//   cont_chir_fit_data_t_epsilon(double aml,double ams,size_t ib,size_t L,dboot_t wfse,dboot_t wofse) : aml(aml),ams(ams),ib(ib),L(L),wfse(wfse),wofse(wofse) {}
// };

// //! perform a fit to the continuum and chiral
// void cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_epsilon> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an);

// class cont_chir_fit_data_t_k
// {
// public:
//   double aml,ams;
//   size_t ib,L;
//   dboot_t y,fse;;
//   cont_chir_fit_data_t_k(double aml,double ams,size_t ib,size_t L,dboot_t y,dboot_t fse) : aml(aml),ams(ams),ib(ib),L(L),y(y),fse(fse) {}
// };

// //! perform a fit to the continuum and chiral
// void cont_chir_fit_k(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_k> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an);

#endif
