#ifndef _HASHEDFUNCTION_HPP
#define _HASHEDFUNCTION_HPP

#include <git_hash.hpp>
#include <raw_file.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

/// Hash a function and provides the interpolation
template <int Degree=3>
struct HashedFunction
{
  /// Path where to save or load
  const string& path;
  
  /// Minimal value of the argument
  const double xMin;
  
  /// Maximal value of the argument
  const double xMax;
  
  /// Number of intervals
  const int nIntervals;
  
  /// Step
  const double dX;
  
  /// Approximation degree
  static constexpr int degree=
    Degree;
  
  static_assert(degree%2,"degree must be odd");
  
  /// Number of interpolation parameters
  static constexpr int nInterpolationPars=
    degree+1;
  
  /// Offset needed to set the spline
  static constexpr int offset=
    (degree-1)/2;
  
  /// Type used to hold the interpolation parameters in each range
  using InterpolationPars=
    array<double,nInterpolationPars>;
  
  /// List of spline parameters in all ranges
  vector<InterpolationPars> interpolationPars;
  
  /// Adds the proper offset to the interval range
  int iIntervalInRange(const int& i) const
  {
    return
      std::max(offset,std::min(i,nIntervals-offset-2));
  }
  
  /// Gets a spline for the given interval
  template <typename F>
  InterpolationPars getSpline(const F& fun,
			      const int& iInterval) const
  {
    /// Central position of the range
    const double x0=
      xMin+iInterval*dX;
    
    vector<double> Al(2*degree+1,0.0);
    InterpolationPars c{};
    
    for(int dI=-(degree-1)/2;dI<=(degree+1)/2;dI++)
      {
	const double xMinusX0=
	  dI*dX;
	
	const double x=
	  xMinusX0+x0;
	
        /// Weight
        double w=
	  1.0;
	
        for(int f=0;f<=2*degree;f++)
          {
            Al[f]+=w;
            if(f<=degree)
	      c[f]+=fun(x)*w;
            w*=xMinusX0;
          }
      }
  
  vector<double> A((degree+1)*(degree+1));
  for(int i=0;i<=degree;i++)
    for(int j=0;j<=degree;j++)
      A[i*(degree+1)+j]=Al[i+j];
  
  //
  
  for(int i=0;i<degree+1;i++)
    {
      double C=A[i*(degree+1)+i];
      for(int j=i;j<degree+1;j++)
	A[i*(degree+1)+j]/=C;
      c[i]/=C;
      
      for(int k=i+1;k<degree+1;k++)
        {
          double C=A[k*(degree+1)+i];
          for(int j=i;j<degree+1;j++)
	    A[k*(degree+1)+j]-=A[i*(degree+1)+j]*C;
          c[k]-=C*c[i];
        }
    }
  
  InterpolationPars res;
  for(int k=degree;k>=0;k--)
    {
      double S=
	0.0;
      
      for(int i=k+1;i<degree+1;i++)
	S+=A[k*(degree+1)+i]*res[i];
      res[k]=c[k]-S;
    }
  
  return
    res;
  }
  
  /// Gets all the spline
  template <typename F>
  vector<InterpolationPars> getSplines(const F& fun) const
  {
    /// Result parameters
    vector<InterpolationPars> res(nIntervals);
    
    /// Check if we can read
    const bool read=
      (path=="")?
      false:
      readFromFile(res);
    
      if(not read)
	{
	  for(int iInterval=0;iInterval<nIntervals;iInterval++)
	    res[iInterval]=
	      getSpline(fun,iIntervalInRange(iInterval));
	  
	  if(path!="")
	    writeToFile(res);
	}
      
    return
      res;
  }
  
  /// Constructor taking the function, range and number of points
  template <typename F>
  HashedFunction(const string& path,
		 const F& fun,
		 const double& xMin,
		 const double& xMax,
		 const int& nPoints) :
    path(path),
    xMin(xMin),
    xMax(xMax),
    nIntervals(nPoints-1),
    dX((xMax-xMin)/nIntervals),
    interpolationPars(getSplines(fun))
  {
  }
  
  /// Evaluate the parameters in the asked point
  double operator()(const double& x) const
  {
    const int iX=
      floor((x-xMin)/dX);
    
    const int i0=
      iIntervalInRange(iX);
    
    const InterpolationPars& coeffs=
      interpolationPars[i0];
    
    double res=
      0.0;
    
    double rx=
      1.0;
    
    const double dist=
      x-(i0*dX+xMin);
    
    for(int i=0;i<=degree;i++)
      {
	res+=coeffs[i]*rx;
	rx*=dist;
      }
    
    return
      res;
  }
  
  /// Repeats the hash
  const char* gitHash=
    GIT_HASH;
  
  /// Terminator
  template <typename F,
	    typename R>
  auto operateOn(const F& f,
		 const R& r) const
  {
    return
      r;
  }
  
  /// Iteratively loop on passed variables
  template <typename F,
	    typename R,
	    typename O,
	    typename...Tail>
  auto operateOn(const F& f,
		 const R& r,
		 const O& o,
		 const char* name,
		 const Tail&...tail) const
  {
    //auto res=
      f(r,o,name);
      
    return
      operateOn(f,r,tail...);
  }
  
  /// Loop on "control vars"
  template <typename F,
	    typename R>
  auto operateOnControlVars(const F& f,
			    const R& r) const
  {
    return
      operateOn(f,
		r,
		gitHash,"GIT_HASH",
		xMin,"xMin",
		xMax,"xMax",
		dX,"dX",
		nIntervals,"nIntervals",
		degree,"Degree");
  }
  
  /// Write to the file
  void writeToFile(const vector<InterpolationPars>& res) const
  {
    /// File where to store
    raw_file_t file(path,"w");
    
    cout<<"Writing to file \""<<path<<"\""<<endl;
    
    operateOnControlVars([&file](const bool& r,
				 const auto& o,
				 const char* name)
    {
      return
	file.bin_write(o);
    },true);
    
    file.bin_write(res.size());
    file.bin_write(res);
  }
  
  /// Reads from the file
  bool readFromFile(vector<InterpolationPars>& res) const
  {
    /// Result of read
    bool status=
      true;
    
    if(not file_exists(path))
      {
	status=false;
	cout<<"File \""<<path<<"\" not found"<<endl;
      }
    
    if(status)
      {
	cout<<"Checking file \""<<path<<"\""<<endl;
	
	raw_file_t file(path,"r");
	
	status=
	  operateOn([&file](const auto& o,
			    const char* name,
			    bool r)
	{
	  if(r)
	    {
		const auto read=
		  file.bin_read<decltype(o)>();
		
		if(read!=o)
		  {
		    r=false;
		    cout<<"Read variable "<<name<<" has value "<<read<<" different form expected "<<o<<endl;
		  }
	      }
	  
	  return
	      r;
	},true);
	
	if(status)
	  {
	    cout<<"File can be read"<<endl;
	    
	    res.resize(file.bin_read<size_t>());
	    file.bin_read(res);
	  }
	else
	  cout<<"File cannot be read"<<endl;
      }
    
    return
      status;
  }
};


#endif
