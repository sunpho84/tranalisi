#ifndef _PERMES_COMBO_HPP
#define _PERMES_COMBO_HPP

#include <REPH/perens.hpp>

//! Holds all info for a combination of quarks
struct permes_combo_t
{
  //! Reference ensemble
  const perens_t& ens;
  
  //! Index of the spectator quark
  const size_t iMs;
  
  //! Index of the forward line quark
  const size_t iMt;
  
  //! Combination relative to the meson combination
  const std::string mesComboTag;
  
  //! Path where to store mesons plots
  std::string mesPlotsPath;
  
  //! Charge of the spectator quark
  const double eS;
  
  //! Charge of the forward line quark
  const size_t eT;
  
  //! Pseudoscalar coupling
  djvec_t ZP;
  
  //! Axial coupling
  djvec_t ZA;
  
  //! Decay constant taken from Pseudoscalar current
  djvec_t fP;
  
  //! Decay constant taken from Axial current
  djvec_t fPbare;
  
  //! Energy
  djvec_t E;
  
  //! Effective mass
  vector<djvec_t> eEff;
  
  //! Pseudoscalar correlation function
  vector<djvec_t> corrPP;
  
  //! T component of axial current
  vector<djvec_t> corrA0P;
  
  //! Z component of axial current
  vector<djvec_t> corrA3P;
  
  //! Load the PP correlation function
  djvec_t load2ptsPP(const size_t iMoms,const size_t iMomt);
  
  //! Load the AP correlation function
  djvec_t load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  //! Load the three points correlation functions
  djvec_t load3pts(const size_t iVA,const size_t iMoms,const size_t iMomt,const size_t iMom0);
  
  //! Computes the axial couplings
  void computeAxialPseudoCouplings();
  
  //! Plot the dispersion relation
  void plotDispRel() const;
  
  //! Load all data
  void load();
  
  //! Constructor
  permes_combo_t(const perens_t& ens,const size_t& iMs,const size_t& iMt,const double& eS,const double& eT) :
    ens(ens),
    iMs(iMs),
    iMt(iMt),
    mesComboTag(combine("Mes_iMs%zu" "_" "iMt%zu",iMs,iMt)),
    eS(eS),
    eT(eT)
  {
    for(auto& q : {&ZP,&ZA,&fP,&fPbare,&E})
      q->resize(ens.nMesKin);
    
    for(auto& q : {&eEff,&corrA0P,&corrA3P,&corrPP})
      q->resize(ens.nMesKin);
    
    load();
  }
  
  //! Perform the 2pts fit
  void fit2pts(const char* fitTag)
  {
    mesPlotsPath=combine("%s/plots/%s/%s",ens.dirPath.c_str(),mesComboTag.c_str(),fitTag);
    mkdir(mesPlotsPath);
    
    computeAxialPseudoCouplings();
    plotDispRel();
  }
};

#endif
