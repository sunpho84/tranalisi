#include <tranalisi.hpp>

//               ~~~~~~~~~~~~~~~~~~~~~~~~~~ T/2
//               }                           |
//     --- -> Q0 --- X --- -> Qt ---         |
//    /                             \        |
//   /                               \       |
// P ---------- <- QS -------------V-A (t)

#include <reph.hpp>

int main(int narg,char **arg)
{
  set_njacks(15);
  
  perens_t ens(".");
  
  const size_t iQuarkS=get<1>(mesons[D_PLUS]);
  const size_t iQuarkT=get<2>(mesons[D_PLUS]);
  
  const size_t iMs=iMassesOfQuarks[iQuarkS][0];
  const size_t iMt=iMassesOfQuarks[iQuarkT][0];
  
  const double eS=quarkCharge[iQuarkS];
  const double eT=quarkCharge[iQuarkT];
  
  permes_combo_t mesCombo(ens,iMs,iMt,eS,eT,"fixed_tfit");
  
  /*
    - loop on all physical mesons
    - | loop on all ensembles
    - | | * read the input
    - | | loop on all meson combo
    - | | | * read 2pts correlators
    - | | |
    - | | | loop on all algorithms to select 2pts fit range
    - | | | | * output to diferent meson combo folder
    - | | | | * do the fit of PP and AP
    - | | | | * check dispersion relation
    - | | | |
    - | | | read V and A correlation functions
    - | | | | loop on all ratios definitions
    - | | | | | build the ratios
    - | | | | | ..
    - | | | | |
    - | | | | | loop on all 3pts fit range
    - | | | | | | fit 3pts
    - | | | |
    - | loop on all systematics so far
    - | | loop on all 8 analysis
    - | | | interpolate to physical meson_corrs
    - | | | fit all ensembles with a multilinear function
   */
  
  return 0;
}
