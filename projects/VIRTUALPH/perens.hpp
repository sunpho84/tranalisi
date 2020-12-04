#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include "nazarioStructures.hpp"

constexpr size_t SMEARED=0;

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
  
  /// Coupling of the local operator
  djack_t ZPloc;
  
  /// Coupling pf the smeared operator
  djack_t ZPsme;
  
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
  
  /// Compute EgT
  double EgT(const size_t iX) const
  {
    return sinh(eG[iX])*(1-exp(-T*eG[iX]));
  }
  
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
    
    // for(size_t icorr=0;icorr<2;icorr++)
    //   for(size_t icomb=0;icomb<(size_t)n.ncomb;icomb++)
    // 	for(size_t alpha=1;alpha<=2;alpha++)
    // 	  {
    // 	    const djvec_t c1=
    // 	      n.threePts[n.index3pts({icorr,icomb,SMEARED,1/*mu*/,alpha,1})];
    // 	    const djvec_t c2=
    // 	      n.threePts[n.index3pts({icorr,icomb,SMEARED,2      ,alpha,1})];
	    
    // 	    const djvec_t c=-(a+b)/sqrt(2);
    // 	    const djvec_t d=(a-b)/sqrt(2);
	    
    // 	    // a.ave_err().write(combine("plots/temp_icorr%zu_icomb%zu_mu%zu_alpha%zu_for_pol.xmg",icorr,icomb,1,alpha));
    // 	    // b.ave_err().write(combine("plots/temp_icorr%zu_icomb%zu_mu%zu_alpha%zu_for_pol.xmg",icorr,icomb,2,alpha));
    // 	    c.ave_err().write(combine("plots/temp_icorr%zu_icomb%zu_alpha%zu_pol1.xmg",icorr,icomb,alpha));
    // 	    d.ave_err().write(combine("plots/temp_icorr%zu_icomb%zu_alpha%zu_pol2.xmg",icorr,icomb,alpha));
    // 	  }
    
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
	    const djvec_t a11=n.get3pts(HA,iComb[iIns],SMEARED,1,1);
	    const djvec_t a22=n.get3pts(HA,iComb[iIns],SMEARED,2,2);
	    
	    c3[iC3({iIns,HA,iX})]=(a11+a22);
	    c3[iC3({iIns,HA,iX})].ave_err().write(combine("plots/temp_a_iX%zu_iIns%zu.xmg",iX,iIns));
	    
	    const djvec_t v12=n.get3pts(HV,iComb[iIns],SMEARED,1,2);
	    const djvec_t v21=n.get3pts(HV,iComb[iIns],SMEARED,2,1);
	    
	    // v12.ave_err().write(combine("plots/temp_v12_icorr1_iX%zu_iIns%zu.xmg",iX,iIns));
	    // v21.ave_err().write(combine("plots/temp_v21_icorr1_iX%zu_iIns%zu.xmg",iX,iIns));
	    
	    c3[iC3({iIns,HV,iX})]=-(v21-v12);
	    c3[iC3({iIns,HV,iX})].ave_err().write(combine("plots/temp_v_iX%zu_iIns%zu.xmg",iX,iIns));
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
    const size_t tmin=14,tmax=27; //hack
    
    djvec_t M(2),Z2(2);
    for(size_t iSmeLev=0;iSmeLev<2;iSmeLev++)
      {
	const djvec_t& c=c2[iSmeLev];
	
	c.ave_err().write(combine("plots/raw_PP_ll_sm%zu.xmg",iSmeLev));
	
	two_pts_fit(Z2[iSmeLev],M[iSmeLev],c,T/2,tmin,tmax,combine("plots/PP_ll_sm%zu.xmg",iSmeLev));
      }
    
    mP=M[1];
    ZPloc=sqrt(Z2[0]);
    ZPsme=Z2[1]/ZPloc;
    
    cout<<"mP: "<<mP.ave_err()<<endl;
    cout<<"Z: "<<ZPloc.ave_err()<<endl;
  }
  
  perens_t() : mu1(0.27903),mu2(0.02363) //hack
  {
    L=32; //hack
    
    read();
    
    fitPP();
    
    norm=djvec_t(T/2+1);
    for(size_t it=0;it<=T/2;it++)
      norm[it]=ZPsme*exp(-mP*it)/(4*mP*mP);
    
    xG.resize(nX);
    
    for(size_t iX=0;iX<nX;iX++)
      xG[iX]=2*eG[iX]/mP;
  }
};

constexpr char perens_t::ingredientsPath[];

#endif
