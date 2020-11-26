#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include "nazarioStructures.hpp"

/// Holds all data for a given ensemble
struct perens_t
{
  static constexpr char ingredientsPath[]="ingredients.dat";
  
  static constexpr size_t nSmeLev=2;
  
  /// Temporal size
  size_t T;
  
  /// Spatial size
  size_t L;
  
  /// Two points pseudoscalar
  vector<djvec_t> c2;
  
  /// Mass
  djack_t mP;
  
  /// Coupling
  djack_t ZP;
  
  /// Three points
  vector<djvec_t> c3;
  
  /// Index to access three point functions
  index_t iC3;
  
  /// Number of different kinetic conditions
  size_t nX;
  
  /// Value of k
  vector<double> kZ;
  
  /// Value of E
  vector<double> eG;
  
  /// Value of X
  djvec_t xG;
  
  /// Normalization for three points
  djvec_t norm;
  
  /// Inserted, forward mass
  const double mu1;
  
  /// Spectator, backward
  const double mu2;
  
  /// Enumberates the 2pts type
  enum{PP,PA0,PA1,PA2,PA3};
  
  static constexpr size_t is=0,il=1; //index to be fetched from inv list
  
  void setIC3()
  {
    iC3.set_ranges({{"Ins",2},{"AV",2},{"X",nX}});
  }
  
  void readIngredients()
  {
    raw_file_t ingredientsFile(ingredientsPath,"r");
    
    ingredientsFile.bin_read(T);
    
    for(size_t iSmeLev=0;iSmeLev<2;iSmeLev++)
      c2[iSmeLev]=djvec_t(T/2+1);
    
    ingredientsFile.bin_read(c2);
    
    ingredientsFile.bin_read(nX);
    
    eG.resize(nX);
    kZ.resize(nX);
    
    ingredientsFile.bin_read(eG);
    ingredientsFile.bin_read(kZ);
    
    setIC3();
    
    c3.resize(iC3.max());
    for(size_t ic3=0;ic3<iC3.max();ic3++)
      c3[ic3]=djvec_t(T/2+1);
      
    ingredientsFile.bin_read(c3);
    
  }
  
  void writeIngredients()
  {
    raw_file_t ingredientsFile(ingredientsPath,"w");
    
    ingredientsFile.bin_write(T);
    
    ingredientsFile.bin_write(c2);
    
    ingredientsFile.bin_write(nX);
    
    ingredientsFile.bin_write(eG);
    
    ingredientsFile.bin_write(kZ);
    
    ingredientsFile.bin_write(c3);
  }
  
  /// Raw read
  void readRaw()
  {
    /// Initialize the reaer
    nazarioReader n(L);
    
    T=n.tmax;
    
    for(size_t iSmeLev=0;iSmeLev<nSmeLev;iSmeLev++)
      c2[iSmeLev]=n.twoPts[iSmeLev].symmetrized();
    
    /// Pairs of combinations which corresponds to a given X
    vector<pair<size_t,size_t>> combOfX;
    
    for(size_t iComb=0;iComb<(size_t)n.ncomb;iComb++)
      {
	const size_t iOppComb=n.iOppComb(iComb);
	
	if(n.comb[iComb].mu1==mu1 and
	   n.comb[iComb].mu2==mu2//  and
	   // iOppComb>iComb
	   )
	  combOfX.push_back({iComb,iOppComb});
      }
    
    nX=combOfX.size();
    cout<<"Number of X: "<<nX<<endl;
    
    // for(size_t iIndepComb=0;iIndepComb<nIndepCombo;iIndepComb++)
    //   cout<<iIndepComb<<" "<<combOfIndepComb[iIndepComb].first<<" "<<combOfIndepComb[iIndepComb].second<<endl;
    
    setIC3();
    
    c3.resize(iC3.max());
    kZ.resize(nX);
    eG.resize(nX);
    
    for(size_t iX=0;iX<nX;iX++)
      {
	const size_t iComb[2]={combOfX[iX].first,combOfX[iX].second};
	
	eG[iX]=n.eG[iComb[0]];
	kZ[iX]=n.kZ[iComb[0]];
	
	for(size_t iIns=0;iIns<2;iIns++)
	  {
	    constexpr size_t SMEARED=1;
	    
	    const djvec_t a11=n.get3pts(HA,iComb[iIns],SMEARED,1,1);
	    const djvec_t a22=n.get3pts(HA,iComb[iIns],SMEARED,2,2);
	    
	    constexpr int coeffA[2]={+1,-1};
	    c3[iC3({iIns,HA,iX})]=(a11+a22)*coeffA[iIns];
	    
	    const djvec_t v12=n.get3pts(HV,iComb[iIns],SMEARED,1,2);
	    const djvec_t v21=n.get3pts(HV,iComb[iIns],SMEARED,2,1);
	    
	    c3[iC3({iIns,HV,iX})]=-(v21-v12);
	  }
      }
  }
  
  void read()
  {
    c2.resize(2);
    
    if(file_exists(ingredientsPath))
      readIngredients();
    else
      {
	readRaw();
	writeIngredients();
      }
  }
  
  /// Fit PP correlation functions
  void fitPP()
  {
    const size_t tmin=15,tmax=28; //hack
    
    djvec_t M(2),Z2(2);
    for(size_t iSmeLev=0;iSmeLev<2;iSmeLev++)
      {
	const djvec_t& c=c2[iSmeLev];
	
	c.ave_err().write(combine("plots/raw_PP_ll_sm%zu.xmg",iSmeLev));
	
	two_pts_fit(Z2[iSmeLev],M[iSmeLev],c,T/2,tmin,tmax,combine("plots/PP_ll_sm%zu.xmg",iSmeLev));
      }
    
    mP=M[1];
    ZP=Z2[1]/sqrt(Z2[0]);
    
    cout<<"mP: "<<mP.ave_err()<<endl;
    cout<<"Z: "<<ZP.ave_err()<<endl;
  }
  
  perens_t() : mu1(0.27903),mu2(0.02363) //hack
  {
    L=32; //hack
    
    read();
    
    fitPP();
    
    norm=djvec_t(T/2+1);
    for(size_t it=0;it<=T/2;it++)
      norm[it]=ZP*exp(-mP*it)/(4*mP*mP);
    
    xG.resize(nX);
    
    for(size_t iX=0;iX<nX;iX++)
      xG[iX]=2*eG[iX]/mP;
  }
};

#endif
