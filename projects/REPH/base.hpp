#ifndef _BASE_HPP
#define _BASE_HPP

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <tranalisi.hpp>

const char VA_tag[2][2]={"V","A"};

const char ST_tag[2]={'S','T'};

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

/////////////////////////////////////////////////////////////////

//! Range in the form [begin;end]
struct Range
{
  size_t begin;
  size_t end;
  
  //Returns the size, including the last point
  size_t size() const
  {
    return end-begin+1;
  }
};

//! Print a range
ostream& operator<<(ostream& os,const Range& r)
{
  return os<<"["<<r.begin<<";"<<r.end<<"]";
}

//! Merges selection if their distance (from their middle) is smaller than their total size
vector<Range> mergeRangesByDistanceSize(vector<Range> ranges,const bool verbose=false)
{
  size_t nMerged;
  size_t iRange;
  
  do
    {
      nMerged=0;
      iRange=0;
      
      while(ranges.size()>1 and iRange<ranges.size()-1)
	{
	  const size_t& begThis=ranges[iRange].begin;
	  const size_t& endThis=ranges[iRange].end;
	  const size_t& begNext=ranges[iRange+1].begin;
	  const size_t& endNext=ranges[iRange+1].end;
	  
	  const size_t nextSize=(endNext-begNext+1);
	  const size_t thisSize=(endThis-begThis+1);
	  
	  const size_t tol=round((nextSize+thisSize+1e-10)/2);
	  const size_t gap=begNext-endThis;
	  
	  if(tol>=gap)
	    {
	      if(verbose) cout<<"Merged "<<ranges[iRange]<<" with "<<ranges[iRange+1]<<", sizes: "<<thisSize<<" "<<nextSize<<", gap: "<<gap<<", tol: "<<tol<<", result";
	      ranges[iRange].end=endNext;
	      ranges.erase(ranges.begin()+iRange+1);
	      if(verbose) cout<<" "<<ranges[iRange]<<endl;
	      
	      nMerged++;
	    }
	  else
	    {
	      if(verbose) cout<<"Not Merged "<<ranges[iRange]<<" with "<<ranges[iRange+1]<<", sizes: "<<thisSize<<" "<<nextSize<<", gap: "<<gap<<", tol: "<<tol<<endl;;
	      iRange++;
	    }
	}
      
      if(verbose)
	{
	  cout<<"NMerged: "<<nMerged<<endl;
	  cout<<endl;
	}
    }
  while(nMerged!=0);
  
  return ranges;
}

//! Find the range in which the quantity is compatible with 0
template <typename TV=djvec_t>
class CompatibilityRangeFinder
{
  //! Reference quantity
  const TV y;
  
  //! Average
  const vec_ave_err_t aveErr;
  
  //! Covariance matrix
  mutable Matrix<double,Dynamic,Dynamic>* _covMatr{nullptr};
  
  //! Considered range beginning
  size_t tSearchMin;
  
  //! Considered range end
  size_t tSearchMax;
  
  //! Significativity of the the quantity
  const vector<double> signif;
  
  //! Selected range
  vector<bool> selected;
  
  //! Plot file
  grace_file_t plot;
  
  //! Count the number of range plots
  size_t rangePlotCounter;
  
  //! Offset between consecutive range plots
  double plotOffset;
  
  //! Verbosity
  bool verbose;
  
public:
  
  //! Return covariance matrix
  const Matrix<double,Dynamic,Dynamic>& covMatr() const
  {
    if(_covMatr==nullptr)
      {
	const size_t& n=y.size();
	
	_covMatr=new Matrix<double,Dynamic,Dynamic>(n,n);
	
	for(size_t i=0;i<n;i++)
	  for(size_t j=0;j<n;j++)
	    (*_covMatr)(i,j)=cov(y[i],y[j]);
      }
    
    return *_covMatr;
  }
  
  //! Correlated chi2 of the fit performed between tMin and tMax
  double ch2CorrFit(const Range& range) const
  {
    const size_t d=range.size();
    if(d==0) CRASH("Null range");
    
    const size_t tMin=range.begin;
    const Matrix<double,Dynamic,Dynamic> cInv=covMatr().block(tMin,tMin,d+1,d+1).inverse();
    
    double Ch2Corr=0;
    for(size_t i=0;i<=d;i++)
      for(size_t j=0;j<=d;j++)
	Ch2Corr+=aveErr[i+tMin].ave()*y[j+tMin].ave()*cInv(i,j);
    
    return Ch2Corr;
  }
  
  //! Check that the range is sensible
  void verifyRange() const
  {
    if(tSearchMax>y.size())
      CRASH("Cannot have %zu as maximal range",tSearchMax);
    
    if(tSearchMax<=tSearchMin)
      CRASH("Cannot have [%zu:%zu) as search range",tSearchMin,tSearchMax);
  }
  
  //! Set the search range to be [frac,1-frac]*size
  CompatibilityRangeFinder& setRangeToConsiderByFraction(const double& frac)
  {
    tSearchMin=round(y.size()*frac);
    tSearchMax=round(y.size()*(1-frac));
    
    verifyRange();
    
    if(verbose) cout<<"Search range: "<<Range{tSearchMin,tSearchMax}<<endl;
    
    return *this;
  }
  
  //! Open the plot path
  CompatibilityRangeFinder& openPlot(const string& path)
  {
    if(plot.is_open())
      plot.close();
    
    rangePlotCounter=0;
    plotOffset=0.0;
    
    if(path!="")
      {
	plot.open(path);
	plot.write_vec_ave_err(aveErr);
      }
    
    return *this;
  }
  
  //! Constructor
  CompatibilityRangeFinder(const TV& y,const string& path="",const bool verbose=false) : y(y),aveErr(y.ave_err()),tSearchMin(1),tSearchMax(y.size()-1),signif(y.significativity()),selected(y.size(),false),verbose(verbose)
  {
    openPlot(path);
  }
  
  //! Set verbosity
  CompatibilityRangeFinder& setVerbose(const bool val)
  {
    verbose=val;
    
    return *this;
  }
  
  //! Destructor
  ~CompatibilityRangeFinder()
  {
    if(_covMatr!=nullptr)
      delete _covMatr;
  }
  
  //! Finds the point which is more compatible with zero, closest to tRef, within nSig (adjusted if needed)
  size_t getClosestCompatiblePointWithinNsigma(const size_t& tRef,double& nSig) const
  {
    size_t tF=0;
    
    bool found=false;
    do
      {
	for(size_t t=tSearchMin;t<tSearchMax;t++)
	  {
	    const bool isComp=signif[t]<nSig;
	    const bool isFirstCheck=tF==0;
	    const bool isCloser=fabs((int)tRef-t)<fabs((int)tF-t);
	    
	    if(isComp and (isFirstCheck or isCloser))
	      {
		found=true;
		tF=t;
	      }
	  }
	
	if(not found)
	  {
	    nSig*=1.1;
	    if(verbose) cout<<"No compatible point found, increasing nSigma to "<<nSig<<endl;
	  }
      }
    while(not found);
    
    if(verbose) cout<<"Most compatible point found: "<<tF<<endl;
    
    return tF;
  }
  
  //! Select the point which is more compatible with zero, closest to tRef, withing nSig (adjusted if needed)
  CompatibilityRangeFinder& selectClosestCompatiblePointWithinNsigma(const size_t& tRef,double& nSig)
  {
    selectPoint(getClosestCompatiblePointWithinNsigma(tRef,nSig));
    
    return *this;
  }
  
  //! Deselect all
  CompatibilityRangeFinder& selectNone()
  {
    fill(selected.begin(),selected.end(),false);
    
    return *this;
  }
  
  //! Select a single point
  CompatibilityRangeFinder& selectPoint(const size_t& t)
  {
    selected[t]=true;
    
    return *this;
  }
  
  //! Find the most compatible point
  size_t getMostCompatiblePoint() const
  {
    size_t tF=tSearchMin;
    
    for(size_t t=tSearchMin+1;t<tSearchMax;t++)
      if(signif[t]<=signif[tF])
	tF=t;
    
    return tF;
  }
  
  //! Select the most compatible point
  CompatibilityRangeFinder& selectMostCompatiblePoint()
  {
    selectPoint(getMostCompatiblePoint());
    
    return *this;
  }
  
  //! Select also all points compatible withing nSig
  CompatibilityRangeFinder& selectAllPointCompatibleWithinNSigma(const double nSig)
  {
    for(size_t t=tSearchMin;t<tSearchMax;t++)
      selected[t]=selected[t] or (signif[t]<=nSig);
    
    return *this;
  }
  
  //! Gets the point with the largest error
  size_t getSelectedPointWithLargestError() const
  {
    if(countSelected()==0)
      CRASH("No point selected");
    
    size_t tF=0;
    
    for(size_t t=tSearchMin;t<tSearchMax;t++)
      if(selected[t])
	if(tF==0 or aveErr[t].err()<aveErr[tF].err())
	  tF=t;
    
    return tF;
  }
  
  //! Gets all points compatible within nSig and errrr not larger of a fraction f
  CompatibilityRangeFinder& selectEnlargingError(const double& f,const double& nSig=1.0)
  {
    const size_t tMaxErr=getSelectedPointWithLargestError();
    if(verbose) cout<<"Point with largest error: "<<tMaxErr<<endl;
    const double eF=f*aveErr[tMaxErr].err();
    
    return
      selectAllPointCompatibleWithinNSigma(nSig).
      trimAllPointsWithErrorLargerThan(eF);
  }
  
  //! Include rightermost point, if compatible within n sigma
  CompatibilityRangeFinder& extendRightWithCompatiblePoints(const double& nSig=1.0)
  {
    size_t t=getSelectionRanges().back().end+1;
    
    while(t<tSearchMax and signif[t]<nSig)
      selected[t++]=true;
    
    return *this;
  }
  
  //! Include lefterrmost point, if compatible within n sigma
  CompatibilityRangeFinder& extendLeftWithCompatiblePoints(const double& nSig=1.0)
  {
    size_t t=getSelectionRanges().front().begin-1;
    
    while(t>tSearchMin and signif[t]<nSig)
      selected[t--]=true;
    
    return *this;
  }
  
  //! Remove from the selection points with too large error
  CompatibilityRangeFinder& trimAllPointsWithErrorLargerThan(const double thresh)
  {
    for(size_t t=tSearchMin;t<tSearchMax;t++)
      if(selected[t] and aveErr[t].err()>thresh)
	selected[t]=false;
    
    return *this;
  }
  
  //! Count selected points
  size_t countSelected() const
  {
    return count(selected.begin()+tSearchMin,selected.begin()+tSearchMax,true);
  }
  
  //! Ranges between fittable ranges
  vector<Range> getSelectionRanges() const
  {
    vector<Range> ranges;
    for(size_t t=tSearchMin;t<=tSearchMax;t++)
      {
	if((t<tSearchMax and selected[t]) and (t==0 or not selected[t-1]))
	  ranges.push_back({t,0});
	if((t>0 and selected[t-1]) and (t>=tSearchMax or not selected[t]))
	  ranges.back().end=t-1;
      }
    
    if(ranges.size()==0)
      {
	cout<<selected<<endl;
      CRASH("No range found!");
      }
    
    //Check that ranges make range
    if(ranges.back().end==0)
      CRASH("Unfinished range");
    
    return ranges;
  }
  
  //! Returns the i-th range
  Range getSelectionRange(const size_t iRange) const
  {
    return getSelectionRanges()[iRange];
  }
  
  //! Select also the passed range
  CompatibilityRangeFinder& selectAlsoRange(const Range& range)
  {
    for(size_t t=range.begin;t<=range.end;t++)
      selected[t]=true;
    
    return *this;
  }
  
  //! Select the passed range
  CompatibilityRangeFinder& selectByRange(const Range& range)
  {
    return
      selectNone().selectAlsoRange(range);
  }
  
  //! Select also the passed ranges
  CompatibilityRangeFinder& selectAlsoRanges(const vector<Range>& ranges)
  {
    for(auto& r : ranges)
      selectAlsoRange(r);
    
    return *this;
  }
  
  //! Select the passed ranges
  CompatibilityRangeFinder& selectByRanges(const vector<Range>& ranges)
  {
    return selectNone().
      selectAlsoRanges(ranges);
  }
  
  //! Merges selection if their distance (from their middle) is smaller than their total size
  CompatibilityRangeFinder& mergeSelectionByDistanceSize()
  {
    return selectByRanges(mergeRangesByDistanceSize(getSelectionRanges(),verbose));
  }
  
  //! Plot the selected points
  CompatibilityRangeFinder& plotSelected()
  {
    const vector<Range> ranges=getSelectionRanges();
    const size_t nRanges=ranges.size();
    
    if(nRanges==0)
      CRASH("No range to be plotted");
    
    if(plotOffset==0.0)
      plotOffset=y[ranges[0].begin].err()*3;
    
    for(size_t iRange=0;iRange<nRanges;iRange++)
      plot.write_line([this](double){return plotOffset/(rangePlotCounter+1);},ranges[iRange].begin-0.25,ranges[iRange].end+0.25,plot.get_line_col_no_increment(),2);
    
    rangePlotCounter++;
    
    return *this;
  }
  
  //! Returns the largest range
  CompatibilityRangeFinder& selectLargestRange()
  {
    const vector<Range> ranges=getSelectionRanges();
    const size_t nRanges=ranges.size();
    size_t largestRange=0;
    
    if(nRanges==0)
      CRASH("No range to be searched");
    
    for(size_t iRange=1;iRange<nRanges;iRange++)
      if(ranges[iRange].size()>ranges[largestRange].size())
	largestRange=iRange;
    
    selectByRange(ranges[largestRange]);
    
    return *this;
  }
  
  //! Returns the compatibility by ch2 test
  vector<double> ch2Compatibility() const
  {
    const size_t& n=y.size();
    vector<double> res(n*n,1e300);
    
    const index_t ind({{"tMin",n},{"tMax",n}});
    
    const size_t D=2;
    for(size_t tMin=tSearchMin;tMin<tSearchMax-D;tMin++)
      for(size_t tMax=tMin+D;tMax<tSearchMax;tMax++)
	res[ind({tMin,tMax})]=ch2CorrFit({tMin,tMax});
    
    return res;
  }
  
  //! Returns the p-value test
  vector<double> getRangePValues() const
  {
    const vector<Range> ranges=getSelectionRanges();
    const size_t nRanges=ranges.size();
    vector<double> pValue(nRanges);
    
    if(nRanges==0)
      CRASH("No range to be searched");
    
    transform(ranges.begin(),ranges.end(),pValue.begin(),[this](const Range& r)
							 {
							   const size_t& n=r.size();
							   const double c=ch2CorrFit(r);
							   return ch2Distr(c/n,n);
							 });
    
    return pValue;
  }
};

// template <bool Verbose=false>
// void compWithZero(size_t& tMin,size_t& tMax,const djvec_t& y,const size_t& tSearchMin,const size_t& tSearchMax,const double incrementDevMin=1.1,const char* plotPath=nullptr)
// {
//   tMin=tMax=0;
  
//   grace_file_t outrange;
//   if(plotPath!=nullptr)
//     outrange.open(plotPath);
  
//   //! Level of compatibility with zero
//   const vector<double> compW0=y.significativity();
  
//   //! Minimal incompatibility
//   const double &compW0Min=*min_element(&compW0[tSearchMin],&compW0[tSearchMax+1]);
//   if(Verbose) cout<<"incomp min: "<<compW0Min<<endl;
  
//   //! Number of stddev for minimal incompatibility
//   const double nDevMin=std::max(compW0Min,incrementDevMin);
//   if(Verbose) cout<<"nDev Minimal: "<<nDevMin<<endl;
  
//   //! Mark where is fittable
//   vector<int> isPlat(compW0.size(),false);
//   std::transform(&compW0[tSearchMin],&compW0[tSearchMax+1],&isPlat[tSearchMin],[nDevMin](double x){return (x<nDevMin);});
//   isPlat[tSearchMin-1]=isPlat[tSearchMax]=false;
  
//   //! Compute error
//   vector<double> err(y.size(),1e300);
//   transform(&y[tSearchMin],&y[tSearchMax+1],&isPlat[tSearchMin],&err[tSearchMin],[](const djack_t& x,const bool& isP){return isP?x.err():1e300;});
  
//   //! take minimal value
//   const auto temp=min_element(&err[tSearchMin],&err[tSearchMax+1]);
//   if(Verbose) cout<<"Min err pos: "<<(temp-&err[0])<<endl;
//   const double minErr=*temp;
  
//   // Remove disgraced points
//   const size_t nQual=2;
//   std::transform(isPlat.begin(),isPlat.end(),err.begin(),isPlat.begin(),[minErr](const bool isP,const double e){return isP and (e<minErr*nQual);});
  
//   //! Margins between fittable ranges
//   vector<size_t> margins;
//   for(size_t t=tSearchMin;t<=tSearchMax;t++)
//     {
//       if(isPlat[t-1] and not isPlat[t])
// 	margins.push_back(t-1);
//       if(isPlat[t] and not isPlat[t-1])
// 	margins.push_back(t);
//     }
  
//   //Check that we have at least a range
//   if(margins.size()==0)
//     CRASH("No compatible point for kinematic");
  
//   //Check that margins make range
//   if(margins.size()%2)
//     CRASH("Size %d should be even",(int)margins.size());
  
//   //Print ranges
//   const size_t nRanges=margins.size()/2;
//   if(Verbose)
//     {
//       cout<<"NRanges: "<<nRanges<<"  ";
//       for(size_t iRange=0;iRange<nRanges;iRange++)
// 	cout<<"["<<margins[iRange*2]<<";"<<margins[iRange*2+1]<<"] ";
//       cout<<endl;
//     }
  
//   if(plotPath)
//     for(size_t iRange=0;iRange<nRanges;iRange++)
//       outrange.write_line([v=y[margins[iRange*2]].err()*iRange/(2*nRanges)](double){return 0.0;},margins[iRange*2]-0.25,margins[iRange*2+1]+0.25);
  
//   size_t iRange;
//   int nMerged;
//   do
//     {
//       nMerged=0;
//       iRange=0;
//       while(margins.size()/2>1 and iRange<margins.size()/2)
// 	{
// 	  const size_t& begThis=margins[iRange+0];
// 	  const size_t& endThis=margins[iRange+1];
// 	  const size_t& begNext=margins[iRange+2];
// 	  const size_t& endNext=margins[iRange+3];
	  
// 	  const size_t nextSize=(endNext-begNext+1);
// 	  const size_t thisSize=(endThis-begThis+1);
	  
// 	  const size_t tol=(nextSize+thisSize)/2;
// 	  const size_t gap=begNext-endThis;
	  
// 	  if(tol>=gap)
// 	    {
// 	      if(Verbose) cout<<"Merged ["<<begThis<<";"<<endThis<<"] with ["<<begNext<<";"<<endNext<<"], sizes: "<<thisSize<<" "<<nextSize<<", gap: "<<gap<<", tol: "<<tol<<", result";
// 	      margins.erase(margins.begin()+iRange+1,margins.begin()+iRange+2+1);
// 	      if(Verbose) cout<<" ["<<begThis<<";"<<endThis<<"]"<<endl;
	      
// 	      nMerged++;
// 	    }
// 	  else
// 	    {
// 	      if(Verbose) cout<<"Not merged ["<<begThis<<";"<<endThis<<"] with ["<<begNext<<";"<<endNext<<"] sizes: "<<thisSize<<" "<<nextSize<<", gap: "<<gap<<", tol: "<<tol<<endl;
// 	      iRange+=2;
// 	    }
// 	}
//       if(Verbose)
// 	{
// 	  cout<<"NMerged: "<<nMerged<<endl;
// 	  cout<<endl;
// 	}
//     }
//   while(nMerged!=0);
  
//   if(plotPath!=nullptr)
//     {
//       for(size_t iRange=0;iRange<margins.size()/2;iRange++)
// 	outrange.write_line([x=y[margins[iRange*2]].err()](double){return x;},margins[iRange*2]-0.25,margins[iRange*2+1]+0.25);
      
//       outrange.write_vec_ave_err(y.ave_err());
//     }
  
//   //select tmin/max taking largest interval
//   for(size_t iRange=0;iRange<margins.size()/2;iRange++)
//     if(tMax-tMin<=margins[iRange*2+1]-margins[iRange*2])
//       {
// 	tMax=margins[iRange*2+1];
// 	tMin=margins[iRange*2];
//       }
  
//   if(Verbose)
//     cout<<" range ["<<tMin<<";"<<tMax<<"]"<<endl;
// }

#endif
