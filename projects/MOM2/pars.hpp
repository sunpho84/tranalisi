#ifndef _PARS_HPP
#define _PARS_HPP

#ifndef EXTERN_PARS
 #define EXTERN_PARS extern
 #define INIT_PARS_TO(...)
#else
 #define INIT_PARS_TO(...) __VA_ARGS__
#endif

#include <tranalisi.hpp>

#include <MOM2/types.hpp>

namespace gaz
{
  enum type_t{PLAQ,TLSYM,IWA};
  const map<string,type_t> decr={{"Plaq",PLAQ},{"tlSym",TLSYM},{"Iwa",IWA}};
  const size_t n=3;
  PROVIDE_DECRYPTER;
}

namespace gf
{
  enum type_t{FEYNMAN,LANDAU};
  const map<string,type_t> decr={{"Feynman",FEYNMAN},{"Landau",LANDAU}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
  const double lambda[2]={1.0,0.0};
}

//! chiral extrapolation
namespace chir_extr
{
  enum type_t{MQUARK,MMES};
  const map<string,type_t> decr={{"MQuark",MQUARK},{"MMes",MMES}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}

//! kind of scheme supported
namespace reno_scheme
{
  enum type_t{RI_MOM,SMOM};
  const map<string,type_t> decr={{"RI_MOM",RI_MOM},{"SMOM",SMOM}};
  const size_t n=2;
  PROVIDE_DECRYPTER;
}

namespace pars
{
  //! allow changing a few pars
  EXTERN_PARS bool can_change_pars INIT_PARS_TO({true});
  //! gauge action
  EXTERN_PARS gaz::type_t act INIT_PARS_TO({gaz::PLAQ});
  //! extrapolation method
  EXTERN_PARS chir_extr::type_t chir_extr_method INIT_PARS_TO({chir_extr::MQUARK});
  //! compute ri Zq
  EXTERN_PARS size_t compute_RI INIT_PARS_TO({0});
  //! compute bilinear vertex
  EXTERN_PARS bool compute_bilinears INIT_PARS_TO({true});
  //! compute mesoleptonic vertex
  EXTERN_PARS bool compute_meslep INIT_PARS_TO({false});
  //! use g2 or g2tilde to correct Oa2
  EXTERN_PARS bool correct_Oa2_using_gtilde INIT_PARS_TO({false});
  //! use p or ptilde to correct Oa2
  EXTERN_PARS bool correct_Oa2_using_ptilde INIT_PARS_TO({false});
  //! Filter for democracy
  EXTERN_PARS double filter_thresh INIT_PARS_TO({1});
  //! assume physical basis propagator
  EXTERN_PARS bool phys_basis INIT_PARS_TO({true});
  //! number of momenta between each print
  EXTERN_PARS size_t print_each_mom INIT_PARS_TO({100});
  //! reno scheme used
  EXTERN_PARS reno_scheme::type_t scheme INIT_PARS_TO({reno_scheme::RI_MOM});
  //! assume twisted run
  EXTERN_PARS bool twisted_run INIT_PARS_TO({true});
  //! list of options to compute deltam
  enum DELTAM_METHOD{FROM_PROP,FROM_CORR};
  //! method to be used to compute deltam
  EXTERN_PARS DELTAM_METHOD deltam_method INIT_PARS_TO({FROM_CORR});
  //! include Zq corrections in bilinears
  EXTERN_PARS bool include_Zq_corrections_in_bilinears INIT_PARS_TO({true});
  //! include self energy and tadpole in bilinears
  EXTERN_PARS bool include_self_energy_in_bilinears INIT_PARS_TO({true});
  //! compute deltam_cr counterterm and subtract
  EXTERN_PARS bool use_deltam_cr_ct INIT_PARS_TO({true});
  //! compute deltam_tm counterterm and subtract
  EXTERN_PARS bool use_deltam_tm_ct INIT_PARS_TO({true});
  //! use QED
  EXTERN_PARS size_t use_QED INIT_PARS_TO({0});
  //! reference p2 to which interpolate, in GeV^2
  EXTERN_PARS double p2ref INIT_PARS_TO({13.0});
  //! range fora2p2 fit
  EXTERN_PARS pair<size_t,size_t> extrapolate_to_0_p2_range INIT_PARS_TO({1,2});
  //! perform the allocate immediately
  EXTERN_PARS bool allocate_immediately INIT_PARS_TO({true});
  //! perform the equiv momenta average immediately
  EXTERN_PARS bool average_equiv_momenta_immediately INIT_PARS_TO({false});
  //! perform the r average immediately
  EXTERN_PARS bool average_r_immediately INIT_PARS_TO({false});
  //! compute mPCAC
  EXTERN_PARS bool report_mPCAC INIT_PARS_TO({false});
  //! subtract the shift induced by QED even if not pole is present
  EXTERN_PARS bool sub_meson_mass_shift_when_no_pole INIT_PARS_TO({true});
  //! store whether we can compute also QED
  EXTERN_PARS bool can_compute_QED INIT_PARS_TO({true});
  //! store whether QED Z are relative or not
  EXTERN_PARS bool Z_QED_are_relative INIT_PARS_TO({true});
  //! plot or not the WI determination of Zs/Zp
  EXTERN_PARS bool plot_Zs_fr_Zp INIT_PARS_TO({true});
  
  //! list of ensembles
  EXTERN_PARS vector<string> ens;
  
  inline void set_use_QED(size_t use)
  {
    if(not can_change_pars) CRASH("Cannot change use_QED");
    use_QED=use;
  }
}

//! return the suffix for QED tag
inline std::string QED_tag_suffix()
{
  return pars::Z_QED_are_relative?"_QED_rel":"_QED";
}

//! set all parameters once forever
void freeze_pars();

#undef INIT_PARS_TO
#undef EXTERN_PARS

#endif
