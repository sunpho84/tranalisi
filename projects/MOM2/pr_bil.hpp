#ifndef _PR_BIL_HPP
#define _PR_BIL_HPP

#include <Dirac.hpp>
#include <meas_vec.hpp>

#include <MOM2/geometry.hpp>
#include <MOM2/prop.hpp>

#ifndef EXTERN_PR_BIL
 #define EXTERN_PR_BIL extern
 #define INIT_PR_BIL_TO(...)
#else
 #define INIT_PR_BIL_TO(...) = __VA_ARGS__
#endif

namespace pr_bil
{
  enum ins{LO,QED, CR_OU,CR_IN, TM_OU,TM_IN, PH_OU,PH_IN, QED_OU,QED_IN, EX};
  EXTERN_PR_BIL vector<size_t>   iins_of_ins;
  EXTERN_PR_BIL vector<ins>      ins_list;
  EXTERN_PR_BIL vector<string>   ins_tag INIT_PR_BIL_TO({"LO","QED","CR_OU","CR_IN","TM_OU","TM_IN","PH_OU","PH_IN","QED_OU","QED_IN","EX"});
  
  //! set all pr_bil insertion
  void set_ins();
  EXTERN_PR_BIL size_t nins;
}

//! bilinears
enum ibil_t{iS,iA,iP,iV,iT};
const string bil_tag="SAPVT";

//! number of bilinears
const size_t nbil=5;

//! holds all ibil_t
const ibil_t ibil_t_list[nbil]={iS,iA,iP,iV,iT};
const vector<vector<size_t>> iG_of_bil={{0},{6,7,8,9},{5},{1,2,3,4},{10,11,12,13,14,15}};

#undef EXTERN_PR_BIL
#undef INIT_PR_BIL_TO

#endif
