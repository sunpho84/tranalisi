#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_JACK
#include <jack.hpp>

void jackknivesFill(const size_t& nConfs,const function<void(const size_t& iConf,const size_t& iClust,const double& weight)>& f)
{
  /// Cluster size
  const double clustSize=(double)nConfs/njacks;
  
  // We loop on the clusters, and inside we loop on the configurations to be included
  for(size_t iClust=0;iClust<njacks;iClust++)
    {
      /// Initial time of the bin
      const double binBegin=iClust*clustSize;
      /// Final time of the bin
      const double binEnd=binBegin+clustSize;
      
      /// Loop time
      double binPos=binBegin;
      do
	{
	  /// Index of the configuration related to the time
	  const size_t iConf=floor(binPos+1e-10);
	  
	  ///Rectangle left point
	  const double beg=binPos;
	  
	  /// Rectangle right point
	  const double end=std::min(binEnd,beg+1.0);
	  
	  /// Rectangle horizontal size
	  const double weight=end-beg;
	  
	  // Perform the operation passing the info
	  f(iConf,iClust,weight);
	  
	  // Updates the position
	  binPos=end;
	}
      while(binEnd-binPos>1e-10);
    }
}
