#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <filesystem>

#include <tranalisi.hpp>

#include "global.hpp"

#include "nazarioStructures.hpp"

size_t useSmeared;
const string compile_time=__TIME__;

/// Holds all data for a given ensemble
struct perens_t
{
  /// Name of the enesemble
  const string name;
  
  /// Returns the string needed for the plots path
  const string plotsPath;
  
  /// Name of ingredients file
  static constexpr char ingredientsPath[]="ingredients.dat";
  
  /// Number of smearing levels
  size_t nSmLev;
  
  /// Temporal size
  size_t T;
  
  /// Spatial size
  size_t L;
  
  /// Two points pseudoscalar
  vector<djvec_t> c2;
  
  /// Mass
  djack_t mP;
  
  /// Coupling of the operator
  vector<djack_t> ZP;
  
  /// Three points
  vector<djvec_t> c3;
  
  /// Index to access three point functions
  index_t iC3;
  
  /// Value of k
  vector<double> kZ;
  
  /// Value of E
  vector<double> eG;
  
  /// Minimal and maximal ranges for PP
  size_t tMinPP,tMaxPP;
  
  /// Minimal and maximal ranges for A
  vector<size_t> tMinA,tMaxA;
  
  /// Minimal and maximal range for V
  vector<size_t> tMinV,tMaxV;
  
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
  double mu1;
  
  /// Spectator, backward
  double mu2;
  
  /// Enumberates the 2pts type
  enum{PP,PA0,PA1,PA2,PA3};
  
  static constexpr size_t is=0,il=1; //index to be fetched from inv list
  
  void setIC3()
  {
    iC3.set_ranges({{"Ins",2},{"AV",2},{"X",nX}});
  }
  
  /// Read ingredients
  void readIngredients()
  {
    const string fullPath=name+"/"+ingredientsPath;
    
    if(std::filesystem::last_write_time(fullPath)<std::filesystem::last_write_time(name+"/data/conf.virtualph.dat"))
      CRASH("ingredients file older than data");
    
    raw_file_t ingredientsFile(fullPath,"r");
    
    const string read_compile_time=ingredientsFile.bin_read<string>();
    if(read_compile_time!=compile_time)
      CRASH("Read compile time \"%s\" not agreeing with current one \"%s\"",read_compile_time.c_str(),compile_time.c_str());
    
    const size_t read_nJacks=ingredientsFile.bin_read<size_t>();
    if(read_nJacks!=njacks)
      CRASH("Read njacks %zu not agreeing with current one %zu",read_nJacks,njacks);
    
    const size_t read_useSmeared=ingredientsFile.bin_read<size_t>();
    if(useSmeared!=read_useSmeared)
      CRASH("Read useSmeared %zu not agreeing with current one %zu",read_useSmeared,useSmeared);
    
    ingredientsFile.bin_read(T);
    ingredientsFile.bin_read(nSmLev);
    cout<<nSmLev<<endl;
    
    c2.resize(nSmLev,djvec_t(T/2+1));
    ZP.resize(nSmLev);
    
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
    raw_file_t ingredientsFile(name+"/"+ingredientsPath,"w");
    
    ingredientsFile.bin_write(compile_time);
    
    ingredientsFile.bin_write(njacks);
    
    ingredientsFile.bin_write(useSmeared);
    
    ingredientsFile.bin_write(T);
    
    ingredientsFile.bin_write(nSmLev);
    
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
    nazarioReader n(name,L);
    
    T=n.tmax;
    nSmLev=n.nqsml;
    ZP.resize(nSmLev);
    c2.resize(nSmLev);
    
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
	    
     	    // a.ave_err().write(combine(plotsPath+"temp_icorr%zu_icomb%zu_mu%zu_alpha%zu_for_pol.xmg",icorr,icomb,1,alpha));
     	    // b.ave_err().write(combine(plotsPath+"temp_icorr%zu_icomb%zu_mu%zu_alpha%zu_for_pol.xmg",icorr,icomb,2,alpha));
    // 	    c.ave_err().write(combine(plotsPath+"temp_icorr%zu_icomb%zu_alpha%zu_pol1.xmg",icorr,icomb,alpha));
    // 	    d.ave_err().write(combine(plotsPath+"temp_icorr%zu_icomb%zu_alpha%zu_pol2.xmg",icorr,icomb,alpha));
    // 	  }
    
    for(size_t iSmeLev=0;iSmeLev<nSmLev;iSmeLev++)
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
    
    /// Number of nX computed
    const size_t computedNx=combOfX.size();
    if(computedNx!=nX)
      CRASH("Expected %zu nX, obtained %zu",nX,computedNx);
    
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
	    const djvec_t a11=n.get3pts(HA,iComb[iIns],useSmeared,1,1
					,combine("%stemp_a11_icorr1_iX%zu_iIns%zu.xmg",plotsPath.c_str(),iX,iIns)
					);
	    const djvec_t a22=n.get3pts(HA,iComb[iIns],useSmeared,2,2
					,combine("%stemp_a22_icorr1_iX%zu_iIns%zu.xmg",plotsPath.c_str(),iX,iIns)
					);
	    
	    // a11.ave_err().write(combine("%stemp_a11_icorr1_iX%zu_iIns%zu.xmg",plotsPath.c_str(),iX,iIns));
	    // a22.ave_err().write(combine("%stemp_a22_icorr1_iX%zu_iIns%zu.xmg",plotsPath.c_str(),iX,iIns));
	    
	    c3[iC3({iIns,HA,iX})]=(a11+a22);
	    // c3[iC3({iIns,HA,iX})].ave_err().write(combine(plotsPath+"temp_a_iX%zu_iIns%zu.xmg",iX,iIns));
	    
	    const djvec_t v12=n.get3pts(HV,iComb[iIns],useSmeared,1,2);
	    const djvec_t v21=n.get3pts(HV,iComb[iIns],useSmeared,2,1);
	    
	    // v12.ave_err().write(combine(plotsPath+"temp_v12_icorr1_iX%zu_iIns%zu.xmg",iX,iIns));
	    // v21.ave_err().write(combine(plotsPath+"temp_v21_icorr1_iX%zu_iIns%zu.xmg",iX,iIns));
	    
	    c3[iC3({iIns,HV,iX})]=-(v21-v12);
	    // c3[iC3({iIns,HV,iX})].ave_err().write(combine(plotsPath+"temp_v_iX%zu_iIns%zu.xmg",iX,iIns));
	  }
      }
  }
  
  /// Read the ingredients or compute from raw
  void readOrComputeIngredients()
  {
    if(file_exists(ingredientsPath))
      readIngredients();
    else
      {
	readRaw();
	writeIngredients();
      }
  }
  
  /// Read all parameters
  void readPars()
  {
    raw_file_t input(name+"/pars.txt","r");
    
    L=input.read<size_t>("L");
    mu1=input.read<double>("mu1");
    mu2=input.read<double>("mu2");
    
    input.expect("TIntsPP");
    tMinPP=input.read<size_t>();
    tMaxPP=input.read<size_t>();
    
    input.expect("TIntsA");
    tMinA.resize(nX);
    tMaxA.resize(nX);
    for(size_t iX=0;iX<nX;iX++)
      {
	input.expect(combine("%zu",iX).c_str());
	tMinA[iX]=input.read<size_t>();
	tMaxA[iX]=input.read<size_t>();
      }
    
    input.expect("TIntsV");
    tMinV.resize(nX);
    tMaxV.resize(nX);
    for(size_t iX=0;iX<nX;iX++)
      {
	input.expect(combine("%zu",iX).c_str());
	tMinV[iX]=input.read<size_t>();
	tMaxV[iX]=input.read<size_t>();
      }
  }
  
  /// Fit PP correlation functions
  void fitPP()
  {
    djvec_t M(nSmLev),Z2(nSmLev);
    for(size_t iSmeLev=0;iSmeLev<nSmLev;iSmeLev++)
      {
	const djvec_t& c=c2[iSmeLev];
	
	// c.ave_err().write(combine(plotsPath+"raw_PP_ll_sm%zu.xmg",iSmeLev));
	
	two_pts_fit(Z2[iSmeLev],M[iSmeLev],c,T/2,tMinPP,tMaxPP,plotsPath+combine("PP_ll_sm%zu.xmg",iSmeLev));
      }
    
    mP=M[0];
    ZP[0]=sqrt(Z2[0]);
    for(size_t iSmLev=0;iSmLev<nSmLev;iSmLev++)
      ZP[iSmLev]=Z2[iSmLev]/ZP[0];
    
    cout<<"mP: "<<mP.ave_err()<<endl;
    cout<<"Z: "<<ZP[0].ave_err()<<endl;
  }
  
  perens_t(const string name) : name(name),plotsPath(name+"/plots/")
  {
    readPars();
    
    readOrComputeIngredients();
    
    fitPP();
    
    norm=djvec_t(T/2+1);
    for(size_t it=0;it<=T/2;it++)
      norm[it]=ZP[useSmeared]*exp(-mP*it)/(4*mP*mP);
    
    xG.resize(nX);
    
    for(size_t iX=0;iX<nX;iX++)
      xG[iX]=2*eG[iX]/mP;
  }
};

constexpr char perens_t::ingredientsPath[];

#endif
