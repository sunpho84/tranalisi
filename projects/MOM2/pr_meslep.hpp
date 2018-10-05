#ifndef _PR_MESLEP_HPP
#define _PR_MESLEP_HPP

#include <Dirac.hpp>

#include <MOM2/prop.hpp>

#ifndef EXTERN_PR_MESLEP
 #define EXTERN_PR_MESLEP extern
 #define INIT_PR_MESLEP_TO(...)
#else
 #define INIT_PR_MESLEP_TO(...) __VA_ARGS__
#endif

namespace meslep
{
  //! holds info on constructing operators
  struct Zop_t
  {
    struct listGl_Gq_t
    {
      const size_t ilistGl; //!< insertions to be taken from the list
      const size_t Gq;      //!< basic gamma on the quark side of the operator...
    };
    
    const vector<listGl_Gq_t> contr;
    const int Qg5_sign;               //!< 1+sign*g5 on the quark side
    const size_t norm;                //!< norm when inserting the operator
    const size_t pnorm;               //!< norm when projecting
  };
  
  EXTERN_PR_MESLEP vector<Zop_t> zops INIT_PR_MESLEP_TO({
      {{{1,1},{2,2},{3,3},{4,4}},-1,1,4},
      {{{1,1},{2,2},{3,3},{4,4}},+1,1,4},
      {{{ 0,0}},-1,1,1},
      {{{ 0,0}},+1,1,1},
      {{{{5,10},{6,11},{7,12},{8,13},{9,14},{10,15}}},+1,2,24}});
  
  const size_t nZop_max=5;
  
  EXTERN_PR_MESLEP size_t nZop INIT_PR_MESLEP_TO(=nZop_max);
  
  inline void set_nZops(const int n)
  {
    if(n>(int)nZop_max or n<0) CRASH("Invalid number of nmeslep_ops, %d",n);
    nZop=n;
  }
  
  EXTERN_PR_MESLEP       vector<size_t>         listGl   INIT_PR_MESLEP_TO(={{ 0, 1, 2, 3, 4,10,11,12,13,14,15}}); //!< list of Gamma to be inserted on lepton side of the operator
  EXTERN_PR_MESLEP       vector<size_t>         &listpGl INIT_PR_MESLEP_TO(=listGl);
  EXTERN_PR_MESLEP       vector<int>            Lg5_sign INIT_PR_MESLEP_TO(={{+1,-1,-1,-1,-1,+1,+1,+1,+1,+1,+1}}); //!< 1+sign*g5 on the quark side
}

namespace pr_meslep
{
  enum ins{LO,QED, NA_OU,NA_IN, CR_OU,CR_IN, TM_OU,TM_IN, PH_OU,PH_IN, QED_OU,QED_IN, EX};
  EXTERN_PR_MESLEP vector<size_t>   iins_of_ins;
  EXTERN_PR_MESLEP vector<ins>      ins_list;
  EXTERN_PR_MESLEP vector<string>   ins_tag INIT_PR_MESLEP_TO({"LO","QED","NA_IN","NA_OU","CR_OU","CR_IN","TM_OU","TM_IN","PH_OU","PH_IN","QED_OU","QED_IN","EX"});
  
  //! set all pr_meslep insertion
  void set_ins();
  EXTERN_PR_MESLEP size_t nins;
}

#undef EXTERN_PR_MESLEP
#undef INIT_PR_MESLEP_TO

#endif
