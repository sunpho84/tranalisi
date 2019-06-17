#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <tranalisi.hpp>

const char VA_tag[2][2]={"V","A"};

//! Holds all info for a given physical quark
using quark_t=tuple<double,vector<size_t>>;

//! List of physical quarks
//!
//! Each entry in the map is connected to a physical quark (Up, Down, etc)
//! Each value is a tuple: the first component is the charge, the second is a vector, listing all masses in the data
map<string,quark_t> quarks;

//! Holds all info for a given meson
using meson_t=pair<string,string>;

//! List of all known mesons
//!
//! Each entry in the map is connected to a physical meson (Pi+, D, etc)
//! Each value is a pair of strings, pointing to the physical quark
map<string,meson_t> mesons;

template <bool Verbose=false>
void compWithZero(size_t& tMin,size_t& tMax,const djvec_t& y,const size_t& tSearchMin,const size_t& tSearchMax,const double incrementDevMin=1.1,const char* plotPath=nullptr)
{
  tMin=tMax=0;
  
  grace_file_t outrange;
  if(plotPath!=nullptr)
    outrange.open(plotPath);
  
  //! Level of compatibility with zero
  const vector<double> compW0=y.significativity();
  cout<<compW0<<endl;
  
  //! Minimal incompatibility
  const double &compW0Min=*min_element(&compW0[tSearchMin],&compW0[tSearchMax+1]);
  if(Verbose) cout<<"incomp min: "<<compW0Min<<endl;
  
  //! Number of stddev for minimal incompatibility
  const double nDevMin=std::max(compW0Min*incrementDevMin,1.0);
  if(Verbose)cout<<"nDev Minimal: "<<nDevMin<<endl;
  
  //! Mark where is fittable
  vector<int> isPlat(compW0.size(),false);
  std::transform(&compW0[tSearchMin],&compW0[tSearchMax+1],&isPlat[tSearchMin],[nDevMin](double x){return (x<nDevMin);});
  isPlat[tSearchMin-1]=isPlat[tSearchMax]=false;
  
  //! Compute error
  vector<double> err(y.size(),1e300);
  transform(&y[tSearchMin],&y[tSearchMax+1],&isPlat[tSearchMin],&err[tSearchMin],[](const djack_t& x,const bool& isP){return isP?x.err():1e300;});
  
  //! take minimal value
  const auto temp=min_element(&err[tSearchMin],&err[tSearchMax+1]);
  cout<<"Min err pos: "<<(temp-&err[0])<<endl;
  const double minErr=*temp;
  
  // Remove disgraced points
  const size_t nQual=2;
  std::transform(isPlat.begin(),isPlat.end(),err.begin(),isPlat.begin(),[minErr](const bool isP,const double e){return isP and (e<minErr*nQual);});
  
  //! Margins bwtween fittable ranges
  vector<size_t> margins;
  for(size_t t=tSearchMin;t<=tSearchMax;t++)
    {
      if(isPlat[t-1] and not isPlat[t])
	margins.push_back(t-1);
      if(isPlat[t] and not isPlat[t-1])
	margins.push_back(t);
    }
  
  //Check that we have at least a range
  if(margins.size()==0)
    CRASH("No compatible point for kinematic");
  
  //Check that margins make range
  if(margins.size()%2)
    CRASH("Size %d should be even",(int)margins.size());
  
  //Print ranges
  const size_t nRanges=margins.size()/2;
  cout<<"NRanges: "<<nRanges<<"  ";
  for(size_t iRange=0;iRange<nRanges;iRange++)
    cout<<"["<<margins[iRange*2]<<";"<<margins[iRange*2+1]<<"] ";
  cout<<endl;
  
  if(plotPath)
    for(size_t iRange=0;iRange<nRanges;iRange++)
      outrange.write_line([v=y[margins[iRange*2]].err()*iRange/(2*nRanges)](double){return 0.0;},margins[iRange*2]-0.25,margins[iRange*2+1]+0.25);
  
  size_t iRange;
  int nMerged;
  do
    {
      nMerged=0;
      iRange=0;
      while(margins.size()/2>1 and iRange<margins.size()/2)
	{
	  const auto begThis=margins.begin()+iRange;
	  const auto endThis=margins.begin()+iRange+1;
	  const auto begNext=margins.begin()+iRange+2;
	  const auto endNext=margins.begin()+iRange+3;
	  
	  const auto nextSize=(endNext-begNext+1);
	  const auto thisSize=(endThis-begThis+1);
	  
	  const size_t tol=(nextSize+thisSize)/2;
	  const size_t gap=begNext-endThis;
	  
	  if(tol>=gap)
	    {
	      if(Verbose) cout<<"Merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"], result";
	      margins.erase(endThis,begNext+1);
	      if(Verbose) cout<<" ["<<*begThis<<";"<<*endThis<<"]"<<endl;
	      
	      nMerged++;
	    }
	  else
	    {
	      if(Verbose) cout<<"Not merged ["<<*begThis<<";"<<*endThis<<"] with ["<<*begNext<<";"<<*endNext<<"] sizes: "<<thisSize<<" "<<nextSize<<", gap: "<<gap<<", tol: "<<tol<<endl;
	      iRange+=2;
	    }
	}
      if(Verbose)
	{
	  cout<<"NMerged: "<<nMerged<<endl;
	  cout<<endl;
	}
    }
  while(nMerged!=0);
  
  if(plotPath!=nullptr)
    {
      for(size_t iRange=0;iRange<margins.size()/2;iRange++)
      outrange.write_line([x=y[margins[iRange*2]].err()](double){return x;},margins[iRange*2]-0.25,margins[iRange*2+1]+0.25);
      
      outrange.write_vec_ave_err(y.ave_err());
    }
  
  //select tmin/max taking largest interval
  for(size_t iRange=0;iRange<margins.size()/2;iRange++)
    if(tMax-tMin<=margins[iRange*2+1]-margins[iRange*2])
      {
	tMax=margins[iRange*2+1];
	tMin=margins[iRange*2];
      }
  
  if(Verbose)
    cout<<" range ["<<tMin<<";"<<tMax<<"]"<<endl;
}

#endif
