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
  // vector<double> naiveCovMatr((T/2+1)*(T/2+1),0.0);
  // vector<double> naiveCorrMatr((T/2+1)*(T/2+1),0.0);
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
	// naiveCorrMatr[s+(T/2+1)*t]=
	//   naiveCorrMatr[t+(T/2+1)*s]=
	//   sxy;
	// naiveCovMatr[s+(T/2+1)*t]=
	//   naiveCovMatr[t+(T/2+1)*s]=
	//   sxy/sqrt(sxx*syy);
      }
  
  for(size_t t=0;t<=T/2;t++)
    for(size_t s=0;s<=T/2;s++)
      // covMatr[s+(T/2+1)*t]=naiveCovMatr[s+(T/2+1)*t]*res[s].err()*res[t].err();
  
  // {
  
  // if constexpr(not preciseCov)
  //   {
  //     for(size_t t=0;t<=T/2;t++)
  // 	for(size_t s=0;s<=T/2;s++)
  // 	  covMatr[t+(T/2+1)*s]=cov(res[t],res[s]);
  //   }
  // else
  //   {

  //   for(size_t nBoots=100;nBoots<100000;nBoots*=2)
  //     {
  //   gen_t<ranlux48> gen(91219392);
  
  // // const size_t nBoots=100000;
  
  // vector<double> ave(T/2+1,0.0);
  // vector<double> c((T/2+1)*(T/2+1));
  // vector<double> b(T/2+1,0.0);
  
  // for(auto& ci : c)
  //   ci=0.0;
  
  // for(size_t iBoot=0;iBoot<nBoots;iBoot++)
  //   {
  //     for(auto& bi : b)
  // 	bi=0.0;
      
  //     for(size_t iConf=0;iConf<nConfs;iConf++)
  // 	{
  // 	  const size_t i=gen.get_int(0,nConfs);
  // 	  for(size_t t=0;t<=T/2;t++)
  // 	    b[t]+=raw[iRaw({i,t})];
  // 	}
      
  //     for(size_t t=0;t<=T/2;t++)
  // 	b[t]/=nConfs;
      
  //     for(size_t t=0;t<=T/2;t++)
  // 	{
  // 	  ave[t]+=b[t];
  // 	  for(size_t s=0;s<=T/2;s++)
  // 	    c[s+(T/2+1)*t]+=b[t]*b[s];
  // 	}
  //   }
  
  // for(size_t t=0;t<=T/2;t++)
  //   {
  //     ave[t]/=nBoots;
  //     for(size_t s=0;s<=T/2;s++)
  // 	c[s+(T/2+1)*t]/=nBoots;
  //   }
  
  // for(size_t t=0;t<=T/2;t++)
  //   for(size_t s=0;s<=T/2;s++)
  //     c[t+(T/2+1)*s]=(c[t+(T/2+1)*s]-ave[s]*ave[t])*nConfs;
  
  // // for(size_t t=0;t<=T/2;t++)
  // //   {
  // //     const double f=sqrt(c[t+(T/2+1)*t]);
  // //     for(size_t s=0;s<=T/2;s++)
  // // 	{
  // // 	  c[t+(T/2+1)*s]/=f;
  // // 	  c[s+(T/2+1)*t]/=f;
  // // 	}
  // //   }
  
  // // for(size_t t=0;t<=T/2;t++)
  // //   for(size_t s=t;s<=T/2;s++)
  // size_t t=1,s=2;
  
  // cout<<"AAA "<<nBoots<<" "<<c[t+(T/2+1)*s]<<" "<<naiveCorrMatr[t+(T/2+1)*s]<<endl;
  // }
  // }
  // exit(1);
  
  return make_tuple(T,res,covMatr);
}

