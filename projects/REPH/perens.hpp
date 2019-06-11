#ifndef _PERENS_HPP
#define _PERENS_HPP

#include <tranalisi.hpp>

#include <REPH/base.hpp>

//! Incapsulate all info and operations
struct perens_t
{
  //! Folder where to find data and create plots
  const std::string dirPath;
  
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
  
  //! Momenta of mesons in all kinematics
  vector<double> pMes;
  
  //! Maximal momentum
  double pMesMax;
  
  //! Read the input file
  void readInput()
  {
    raw_file_t fin(combine("%s/jacks/input.txt",dirPath.c_str()),"r");
    
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
    
    pMes.resize(nMesKin);
    
    pMesMax=0;
    for(size_t iKin=0;iKin<nMesKin;iKin++)
      {
	const vector<size_t> c=indMesKin(iKin);
	
	pMes[iKin]=2*M_PI*(moms[c[1]][2]-moms[c[0]][2])/L;
	
	pMesMax=std::max(pMesMax,fabs(pMes[iKin]));
      }
  }
  
  //! Constructor
  perens_t(const std::string dirPath) : dirPath(dirPath)
  {
    readInput();
    
    setKinematics();
  }
};

#endif
