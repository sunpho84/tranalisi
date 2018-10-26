#ifndef _SIGMA_HPP
#define _SIGMA_HPP

#include <tranalisi.hpp>

#ifndef EXTERN_SIGMA
 #define EXTERN_SIGMA extern
 #define INIT_SIGMA_TO(...)
#else
 #define INIT_SIGMA_TO(...) __VA_ARGS__
#endif

namespace sigma
{
  enum ins{LO,  CR,  TM,  PH,  QED,  RI,  RI_QED,  RI_VT,  RI_VX,  RI_VY,  RI_VZ};
  EXTERN_SIGMA vector<ins>     ins_list;
  EXTERN_SIGMA vector<size_t>  iins_of_ins;
  EXTERN_SIGMA vector<string>  ins_tag INIT_SIGMA_TO({"LO","CR","TM","PH","QED","RI","RI_QED","RI_VT","RI_VX","RI_VY","RI_VZ"});
  
  //! set all sigma insertions
  void set_ins();
  EXTERN_SIGMA size_t nins;
  
  enum proj{SIGMA1,SIGMA2,SIGMA3,
	    SIGMA_RI,
	    SIGMA_RI_VT,SIGMA_RI_VX,SIGMA_RI_VY,SIGMA_RI_VZ};
  EXTERN_SIGMA vector<proj>    proj_list;
  EXTERN_SIGMA vector<size_t>  iproj_of_proj;
  EXTERN_SIGMA vector<string>  proj_tag INIT_SIGMA_TO({"sigma1","sigma2","sigma3",
	"sigma_RI",
	"sigma_RI_VT","sigma_RI_VX","sigma_RI_VY","sigma_RI_VZ"});
  
  //! set all sigma proj
  void set_proj();
  EXTERN_SIGMA size_t nproj;
}

#undef EXTERN_SIGMA
#undef INIT_SIGMA_TO

#endif
