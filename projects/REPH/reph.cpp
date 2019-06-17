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
  
  for(size_t iMes=0;iMes<nMesons;iMes++)
    {
      //! Meson to study
      const meson_t& mes=mesons[iMes];
      
      cout<<"Meson: "<<get<0>(mes)<<endl;
      
      const size_t iQuarkS=get<1>(mes);
      const size_t iQuarkT=get<2>(mes);
      
      //! Charge of the spectator and forward line
      const double eS=quarkCharge[iQuarkS];
      const double eT=quarkCharge[iQuarkT];
      
      //! Loop on all meson combos
      const size_t nMs=iMassesOfQuarks[iQuarkS].size();
      const size_t nMt=iMassesOfQuarks[iQuarkT].size();
      const index_t indMesCombo({{"iMs",nMs},{"iMt",nMt}});
      
      //! Holds the 2pts info for each meson combination
      vector<permes_combo_t> mesCombos;
      const size_t nMesCombos=indMesCombo.max();
      
      //! Setup all meson combos
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	{
	  const vector<size_t> c=indMesCombo(iMesCombo);
	  const size_t iMs=iMassesOfQuarks[iQuarkS][c[0]];
	  const size_t iMt=iMassesOfQuarks[iQuarkT][c[1]];
	  mesCombos.emplace_back(ens,iMs,iMt,eS,eT);
	}
      
      for(size_t iMesCombo=0;iMesCombo<nMesCombos;iMesCombo++)
	{
	  mesCombos[iMesCombo]
	    .fit2pts("fixed_tfit")
	    .prepare3ptsNormalization()
	    .chooseTint()
	    .fit3pts("tfixed_tfit")
	    .plotFf();
	}
    }
  
  /*
    - * loop on all physical mesons
    - | loop on all ensembles
    - | | * read the input
    - | | * loop on all meson combo
    - | | | * read 2pts correlators
    - | | |
    - | | | loop on all algorithms to select 2pts fit range
    - | | | | * output to diferent meson combo folder
    - | | | | * do the fit of PP and AP
    - | | | | * check dispersion relation
    - | | | |
    - | | | * read V and A correlation functions
    - | | | | loop on all ratios definitions
    - | | | | | * build the ratios
    - | | | | | ..
    - | | | | |
    - | | | | | loop on all 3pts fit range
    - | | | | | | *fit 3pts
    - | | | |
    - | loop on all systematics so far
    - | | loop on all 8 analysis
    - | | | interpolate to physical meson_corrs
    - | | | fit all ensembles with a multilinear function
   */
  
  return 0;
}
