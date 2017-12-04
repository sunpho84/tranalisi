#ifndef _INGREDIENTS_HPP
#define _INGREDIENTS_HPP

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <read.hpp>
#include <timings.hpp>

#ifndef EXTERN_INGREDIENTS
 #define EXTERN_INGREDIENTS extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) = __VA_ARGS__
#endif

EXTERN_INGREDIENTS index_t imom_ind; //!< index of imom (trivial)
EXTERN_INGREDIENTS index_t r_imom_ind; //!< index of r,imom combo
EXTERN_INGREDIENTS index_t im_r_imom_ind; //!< index of im,r,imom combo
EXTERN_INGREDIENTS index_t indep_imom_ind; //!< index of indep imom (trivial)
EXTERN_INGREDIENTS index_t r_indep_imom_ind; //!< index of r,indep imom combo
EXTERN_INGREDIENTS index_t im_r_indep_imom_ind; //!< index of im,r,indep imom combo
EXTERN_INGREDIENTS index_t im_r_im_r_igam_ind;
EXTERN_INGREDIENTS index_t r_r_iZbil_ind;
EXTERN_INGREDIENTS index_t im_r_im_r_iZbil_ind;
EXTERN_INGREDIENTS index_t iZbil_imom_ind;
EXTERN_INGREDIENTS index_t im_r_im_r_iZbil_imom_ind;
EXTERN_INGREDIENTS index_t im_r_ijack_ind;
EXTERN_INGREDIENTS index_t im_r_ijackp1_ind;

struct ingredients_t
{
  djvec_t deltam_cr;
  
  //Zq for all moms, with and without EM
  djvec_t Zq_allmoms;
  djvec_t Zq_sig1_allmoms;
  djvec_t Zq_sig1_EM_allmoms;
  
  //projected bilinears with and without EM
  djvec_t pr_bil_mom;
  djvec_t pr_bil_QED_mom;
  
  //! Constructor
  ingredients_t();
  
  //! compute deltam_cr
  void compute_deltam_cr();
  
  //! fill all the ingredients according to ri-mom
  void ri_mom();
  
  //! read from binary
  void bin_read(raw_file_t &file);
  
  //! write to binary
  void bin_write(raw_file_t &file) const;
  
  //! read from a path
  inline void bin_read(const string &path)
  {
    raw_file_t fin(path,"r");
    this->bin_read(fin);
  }
  
  //! write to a path
  inline void bin_write(const string &path) const
  {
    raw_file_t fout(path,"w");
    this->bin_write(fout);
  }
};

#undef EXTERN_INGREDIENTS
#undef INIT_TO

#endif
