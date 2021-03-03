#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include <REPH/base.hpp>

#include <regex>

//! Incapsulate all info and operations
struct perens_t
{
  //! Folder where to find data and create plots
  const std::string dirPath;
  
  //! Index of beta
  size_t iBeta;
  
  //! Index of the ultimate analysis
  size_t iUlt;
  
  //! Spatial size
  int L;
  
  //! Temporal size
  int T;
  
  //! Total spatial volume
  int spatVol;
  
  //! Number of masses
  size_t nMass;
  
  //! List of mass
  vector<double> mass;
  
  //! Number of momenta
  size_t nMoms;
  
  //! List of momenta
  vector<array<double,3>> moms;
  
  //! Index spanning all momenta combination
  index_t indMesKin;
  
  //! Number of momenta combination
  size_t nMesKin;
  
  //! Index spanning all decay momenta combination
  index_t indDecKin;
  
  //! Returnt he tag of decKin ind
  string decKinTag(const size_t iMoms,const size_t iMomt,const size_t iMom0) const
  {
    return combine("iMomS%zu_iMomT%zu_iMomO%zu",iMoms,iMomt,iMom0);
  }
  
  //! Wrapper to the function returning the tag
  string decKinTag(const size_t iDecKin) const
  {
    const vector<size_t> c=indDecKin(iDecKin);
    
    return decKinTag(c[0],c[1],c[2]);
  }
  
  //! Number of decay momenta combination
  size_t nDecKin;
  
  //! Momenta of mesons in all kinematics
  vector<double> pMes;
  
  //! Maximal momentum
  double pMesMax;
  
  //! State whether or not to consider to study the decay
  vector<bool> considerDec;
  
  //! State whether there is a symmetric decay
  vector<bool> hasSymmDec;
  
  //! Symmetric of the decay
  vector<int> symmOfDec;
  
  //! Smallest transferred momentum
  vector<int> iDecSmallestKin;
  
  //! Index of meson kinematic of the decay
  vector<int> iMesKinOfDecKin;
  
  //! Energy of the photon
  vector<double> Eg;
  
  //! Finite size correction for the photon energy
  double EgT(const size_t iDecKin) const
  {
    return sinh(Eg[iDecKin])*(1-exp(-T*Eg[iDecKin]));
  }
  
  //! Load a correlation keeping track of it
  djvec_t loadCorr(const string& path,const size_t i,const char* tag=nullptr) const
  {
    if(tag!=nullptr)
      {
	const string loadedListPath=dirPath+"/loaded"+tag+"pts.txt";
	
	static map<string,raw_file_t> m;
	if(m.find(loadedListPath)==m.end())
	  m.emplace(make_pair(loadedListPath,raw_file_t{loadedListPath,"w"}));
	
	m[loadedListPath].printf("%s %zu\n",path.c_str(),i);
      }
    
    return read_djvec(path,T,i);
  }
  
  //! Find X given a mass and kinematic
  double getX(const double m,const double thS,const double thT,const double th0) const
  {
    const double pi=2*M_PI*(th0-thS)/L;
    const double en=latt_en_1D(m,pi);
    const double kHat=2*M_PI*(th0-thT)/L;
    const double eg=2*asinh(fabs(kHat)/2);
    const double pk=en*eg-pi*kHat;
    const double x=2*pk/sqr(en);
    
    return x;
  }
  
  //! Momentum of the decaying pair of lepton
  vector<double> kDec;
  
  //! Lattice version of the decaying momentum
  vector<double> kHatDec;
  
  //! Turns a char into index of beta
  static size_t decryptBeta(const char& b)
  {
    //! Translation from ABD to ib
    const map<char,size_t> betaMap({{'A',0},{'B',1},{'D',2}});
    
    //! Search the character
    auto f=betaMap.find(b);
    if(f==betaMap.end()) CRASH("Unknown %c",b);
    
    return f->second;
  }
  
  //! Gets the "ultimate" index
  static size_t decryptUlt(const string& name)
  {
    //! List of known ultimate input
    static array<string,15> ultMap{"A30.32","A40.32","A50.32","A40.24","A60.24","A80.24","A100.24","B25.32","B35.32","B55.32","B75.32","B85.24","D15.48","D20.48","D30.48"};
    
    //! Search the position of name
    const auto pos=find_if(ultMap.begin(),ultMap.end(),[&name](const string& ult)
						       {
							 return regex_match(name,regex("(.*)("+ult+")([^]*)") );
						       });
    
    if(pos!=ultMap.end())
      return distance(ultMap.begin(),pos);
    else
      return 0;
  }
  
  //! Read the input file
  void readInput()
  {
    raw_file_t fin(combine("%s/jacks/input.txt",dirPath.c_str()),"r");
    
    const string absPath=absolute_path(dirPath);
    cout<<"Absolute path: "<<absPath<<endl;
    
    const string basePath=basename(absPath);
    iBeta=decryptBeta(basePath[0]);
    
    iUlt=decryptUlt(basePath);
    cout<<"Ultimate analysis position: "<<iUlt<<endl;
    
    L=fin.read<size_t>("L");
    spatVol=L*L*L;
    T=fin.read<size_t>("T");
    
    nMass=fin.read<size_t>("NMass");
    mass.resize(nMass);
    for(size_t iMass=0;iMass<nMass;iMass++)
      mass[iMass]=fin.read<double>();
    
    nMoms=fin.read<size_t>("NMoms");
    moms.resize(nMoms);
    for(size_t iMom=0;iMom<nMoms;iMom++)
      for(size_t mu=0;mu<3;mu++)
	moms[iMom][mu]=fin.read<double>();
  }
  
  //! Set all kinematics
  void setKinematics()
  {
    indMesKin.set_ranges({{"mom1",nMoms},{"mom2",nMoms}});
    nMesKin=indMesKin.max();
    
    indDecKin.set_ranges({{"iMoms",nMoms},{"iMomt",nMoms},{"iMom0",nMoms}});
    nDecKin=indDecKin.max();
    
    pMes.resize(nMesKin);
    
    pMesMax=0;
    for(size_t iKin=0;iKin<nMesKin;iKin++)
      {
	const vector<size_t> c=indMesKin(iKin);
	
	pMes[iKin]=2*M_PI*(moms[c[1]][2]-moms[c[0]][2])/L;
	
	pMesMax=std::max(pMesMax,fabs(pMes[iKin]));
	
      }
    
    resizeListOfContainers({&considerDec,&hasSymmDec},indDecKin.max());
    resizeListOfContainers({&symmOfDec,&iMesKinOfDecKin,&iDecSmallestKin},indDecKin.max(),-1);
    resizeListOfContainers({&Eg,&kDec,&kHatDec},indDecKin.max());
    
    ofstream kinematics(dirPath+"/kinematics.txt");
    ofstream symmetrics(dirPath+"/symmetrics.txt");
    kinematics.precision(6);
    kinematics<<std::fixed;
    kinematics<<"Moms\tMomt\tMom0\tk\t\tp"<<endl;
    
    for(size_t iDecKin=0;iDecKin<indDecKin.max();iDecKin++)
      {
     	const vector<size_t> cDec=indDecKin(iDecKin);
	const size_t iMoms=cDec[0],iMomt=cDec[1],iMom0=cDec[2];
	considerDec[iDecKin]=(iMomt!=iMom0);
	
	vector<size_t> c0=cDec;
	c0[1]=c0[2];
	iDecSmallestKin[iDecKin]=indDecKin(c0);
	
	const size_t possSymm=indDecKin({iMomt,iMoms,iMomt});
	hasSymmDec[iDecKin]=(iMoms==iMom0) and (possSymm!=iDecKin);
	if(hasSymmDec[iDecKin])
	  symmOfDec[iDecKin]=possSymm;
	
	if((not hasSymmDec[iDecKin]) or (AVERAGE_SYMMETRIC and (size_t)symmOfDec[iDecKin]<iDecKin))
	  {
	    if(ONLY_SYMMETRIC)
	      considerDec[iDecKin]=false;
	  }
	else
	  symmetrics<<iDecKin<<" "<<symmOfDec[iDecKin]<<endl;
	
	kDec[iDecKin]=2*M_PI*(moms[iMom0][2]-moms[iMomt][2])/L;
	kHatDec[iDecKin]=2*sin(kDec[iDecKin]/2);
	
    	Eg[iDecKin]=2*asinh(fabs(kHatDec[iDecKin])/2);
	iMesKinOfDecKin[iDecKin]=indMesKin({iMoms,iMom0});
	
	const vector<size_t> c=indMesKin(iMesKinOfDecKin[iDecKin]);
	double pMes=2*M_PI*(moms[c[1]][2]-moms[c[0]][2])/L;
	
	
    const double en=latt_en_1D(0.94,pMes);
    const double pk=en*Eg[iDecKin]-pMes*kDec[iDecKin];
    const double x=2*pk/sqr(en);

	kinematics<<iMoms<<"\t"<<iMomt<<"\t"<<iMom0<<"\t"<<kDec[iDecKin]<<"\t"<<pMes<<"\t"<<x<<endl;
      }
  }
  
  //! Constructor
  perens_t(const std::string dirPath) : dirPath(dirPath)
  {
    readInput();
    
    setKinematics();
  }
};

//! Reads the list of ensembles to be used
vector<perens_t> readEnsList()
{
  //! List returned
  vector<perens_t> output;
  
  //! Input file
  raw_file_t input("ensemble_list.txt","r");
  
  //! Total number of ensembles
  const size_t nEns=input.read<size_t>("NEns");
  output.reserve(nEns);
  
  for(size_t iEns=0;iEns<nEns;iEns++)
    {
      const string name=input.read<string>();
      if(ensMap.find(name)==ensMap.end())
	CRASH("Unable to find ensemble %s",name.c_str());
      output.emplace_back(name);
    }
    
  return output;
}

//! Gets the list of quark mass for each quark
map<string,vector<double>> getMassList(const perens_t& e)
{
  //! Result
  map<string,vector<double>> am;
  
  for(auto& q : quarkList)
    for(auto& i : get<2>(q.second))
      am[q.first].push_back(e.mass[i]);
  
  return am;
}

#endif
