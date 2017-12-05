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

EXTERN_INGREDIENTS index_t im_r_im_r_igam_ind;
EXTERN_INGREDIENTS index_t r_r_iZbil_ind;
EXTERN_INGREDIENTS index_t im_r_ijack_ind;
EXTERN_INGREDIENTS index_t im_r_ijackp1_ind;

//! kind of scheme supported
enum scheme_t{RI_MOM};

//! reno scheme used
EXTERN_INGREDIENTS scheme_t scheme;

//! passed to force recomputing Zbil
[[maybe_unused]] static const bool FORCE_RECOMPUTE_ZBIL=true;

struct ingredients_t
{
  size_t _nm;
  size_t _nr;
  vector<double> _am;
  
  vector<vector<size_t>> mom_combo; //!< list of momenta, relatively to glb list
  
  //Zq, with and without EM
  djvec_t Zq;
  djvec_t Zq_sig1;
  djvec_t Zq_sig1_EM;
  
  //projected bilinears with and without EM
  djvec_t pr_bil;
  djvec_t pr_bil_QED;
  
  //bilinear Z
  bool Zbil_computed{false};
  djvec_t Zbil;
  djvec_t Zbil_QED;
  
  index_t imom_ind;
  index_t r_imom_ind;
  index_t im_r_imom_ind;
  index_t indep_imom_ind;
  index_t r_indep_imom_ind;
  index_t im_r_indep_imom_ind;
  index_t iZbil_imom_ind;
  index_t im_r_ind;
  index_t im_r_im_r_iZbil_ind;
  index_t im_r_im_r_iZbil_imom_ind;
  
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
  
  //! create using files
  void create_from_scratch(const string ingredients_path="ingredients.dat");
  
  //! set momenta for ri-moms
  void set_ri_mom_moms();
  
  //! set basic pars
  void set_pars_for_scratch();
  
  //! set all indices
  void set_indices();
  
  //! computes all bilinear Z
  void compute_Zbil();
  
  //allocate all vectors
  void allocate();
  
  //! extrapolate to chiral limit
  ingredients_t chir_extrap() const;
  
  //! subtract O(a^2)
  ingredients_t subtract_Oa2(const bool recompute_Zbil=false) const;
  
  //! evolve to a common scale
  ingredients_t evolve() const;
  
  //! average equivalent momenta
  ingredients_t average_equiv_momenta(const bool recompute_Zbil=false) const;
};

#undef EXTERN_INGREDIENTS
#undef INIT_TO

#endif
