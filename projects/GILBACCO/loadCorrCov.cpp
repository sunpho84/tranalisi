#include <tranalisi.hpp>

/// Loads a correlator and compute covariance
///
/// We don't keep in count the correlation of the data, so the
/// covariance matrix as generated from the boostrap is
/// underestimated. So we pass to the correlation matrix, and then we
/// renormalize with the diagonal part taken from jackknive estimate
/// of the error
auto loadCorrCov(const string& path,
		 const bool computeCov)
{
  /// File needed to read data
  raw_file_t fin(path,"r");
  
  /// Gets the number of confs
  const size_t nConfs=fin.bin_read<size_t>();
  
  /// Read and drops the number of hits
  [[maybe_unused]] const size_t nHits=
		     fin.bin_read<size_t>();
  
  // Read the time extention
  const size_t T=fin.bin_read<size_t>();
  
  /// Index to access the raw data
  index_t iRaw({{"nConfs",nConfs},{"T/2+1",T/2+1}});
  
  /// Raw data
  vector<double> raw(iRaw.max());
  fin.bin_read(raw);
  
  /// Performs the jackknife resampling of the correlator
  djvec_t res(T/2+1);
  jackknivesFill(nConfs,
		 [&](const size_t& iConf,
		     const size_t& iJack,
		     const double& w)
		 {
		   for(size_t t=0;t<=T/2;t++)
		     res[t][iJack]+=w*raw[iRaw({iConf,t})];
		 });
  res.clusterize((double)nConfs/njacks);
  
  /// Covariance matrix
  vector<double> covMatr((T/2+1)*(T/2+1),0.0);
  
  for(size_t t=0;t<=T/2;t++)
    for(size_t s=t;s<=T/2;s++)
      {
	double sx=0,sy=0,sxy=0,sxx=0,syy=0;
	for(size_t iConf=0;iConf<nConfs;iConf++)
	  {
	    const double& x=raw[iRaw({iConf,t})];
	    const double& y=raw[iRaw({iConf,s})];
	    
	    sx+=x;
	    sy+=y;
	    sxx+=x*x;
	    syy+=y*y;
	    sxy+=x*y;
	  }
	for(double* i : {&sx,&sy,&sxx,&syy,&sxy})
	  (*i)/=nConfs;
	  
	sxx-=sx*sx;
	syy-=sy*sy;
	sxy-=sx*sy;
	
	covMatr[s+(T/2+1)*t]=
	  covMatr[t+(T/2+1)*s]=
	  sxy/sqrt(sxx*syy)*res[s].err()*res[t].err();
      }
  
  return make_tuple(T,res,covMatr);
}

