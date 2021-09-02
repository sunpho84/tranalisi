#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_JACK
#include <jack.hpp>

void jackknivesFill(const size_t& nConfs,const function<void(const size_t& iConf,const size_t& iClust,const double& weight)>& f)
{
  const double clustSize=(double)nConfs/njacks;
  
  for(size_t iClust=0;iClust<njacks;iClust++)
    {
      const double confBegin=iClust*clustSize;
      const double confEnd=confBegin+clustSize;
      
      double curConf=confBegin;
      do
	{
	  const size_t iConf=floor(curConf);
	  const double beg=curConf;
	  const double end=std::min(confEnd,beg+1.0);
	  const double weight=end-beg;
	  
	  f(iConf,iClust,weight);
	  
	  curConf=end;
	}
      while(confEnd-curConf>1e-10);
    }
}
