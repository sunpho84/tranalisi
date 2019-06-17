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
  
  //! Path where to store decay plots
  std::string decPlotsPath;
  
  //! Charge of the spectator quark
  const double eS;
  
  //! Charge of the forward line quark
  const double eT;
  
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
  
  //! Difference of energies between initial and final state
  djvec_t dEdec;
  
  //! Quadrimomentum product
  djvec_t PKdec;
  
  //! Form factor independent variable
  djvec_t X;
  
  //! Normalization to be used for 3pts
  vector<djvec_t> normaliz;
  
  //! Pseudoscalar correlation function
  vector<djvec_t> corrPP;
  
  //! T component of axial current
  vector<djvec_t> corrA0P;
  
  //! Z component of axial current
  vector<djvec_t> corrA3P;
  
  //! Decay correlators for V and A
  vector<djvec_t> corrPX[2];
  
  //! Decay correlators for V
  vector<djvec_t>& corrPXV=corrPX[0];
  
  //! Decay correlators for A
  vector<djvec_t>& corrPXA=corrPX[1];
  
  //! Time interval for 3pts fit
  vector<std::pair<size_t,size_t>> tint3pts[2];
  
  //! Form factors for V and A
  djvec_t ff[2];
  
  //! Load the PP correlation function
  djvec_t load2ptsPP(const size_t iMoms,const size_t iMomt);
  
  //! Load the AP correlation function
  djvec_t load2ptsAP(const size_t iMoms,const size_t iMomt,const size_t iGamma);
  
  //! Load the three points correlation functions
  djvec_t load3pts(const size_t iVA,const size_t iMs,const size_t iMt,const size_t iMoms,const size_t iMomt,const size_t iMom0);
  
  //! Computes the axial couplings
  void computeAxialPseudoCouplings();
  
  //! Plot the dispersion relation
  void plotDispRel() const;
  
  //! Load all 2pts
  void load2pts();
  
  //! Load all 3pts
  void load3pts();
  
  //! Load all data
  void load()
  {
    load2pts();
    load3pts();
  }
  
  //! Constructor
  permes_combo_t(const perens_t& ens,const size_t& iMs,const size_t& iMt,const double& eS,const double& eT) :
    ens(ens),
    iMs(iMs),
    iMt(iMt),
    mesComboTag(combine("Mes_iMs%zu" "_" "iMt%zu",iMs,iMt)),
    eS(eS),
    eT(eT)
  {
    resizeListOfContainers({&ZP,&ZA,&fP,&fPbare,&E},ens.nMesKin);
    
    resizeListOfContainers({&eEff,&corrA0P,&corrA3P,&corrPP},ens.nMesKin);
    
    resizeListOfContainers({&ff[0],&ff[1]},ens.nDecKin);
    
    resizeListOfContainers({&tint3pts[0],&tint3pts[1]},ens.nDecKin,std::pair<size_t,size_t>{0,0});
    
    resizeListOfContainers({&corrPXA,&corrPXV},ens.nDecKin);
    
    resizeListOfContainers({&dEdec,&PKdec,&X},ens.nDecKin);
    
    normaliz.resize(ens.indDecKin.max(),djvec_t{(size_t)ens.T/2+1});
    
    load();
  }
  
  //! Prepare the 3pts normalization
  permes_combo_t& prepare3ptsNormalization();
  
  //! Perform the 2pts fit
  permes_combo_t& fit2pts(const char* fitTag);
  
  permes_combo_t& chooseTint();
  
  //! Perform the 3pts fit
  permes_combo_t& fit3pts(const char* fitTag);
  
  //! Plot ff
  permes_combo_t& plotFf();
};

#endif
