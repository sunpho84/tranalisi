#include "../PIPI/common.hpp"

int main()
{
  set_njacks(50);
  
  const std::string dataPath="out";
  const std::vector<std::string> confs=getConfs("confsDataList.dat",dataPath,"");
  const size_t nConfs=confs.size();
  cout<<"NConfs: "<<nConfs<<endl;
  
  const auto rawData=
    getRaw("RUN0.dat",
	   "mes_contr_HH",
	   {""},
	   T,
	   dataPath,
	   confs,
	   {"P5P5","V0P5","V1V1"});
  
  const index_t idx({{"t",T},{"conf",nConfs}});
  
  auto getRaw=
    [&idx,
     &nConfs](const auto& rawData,
	      const std::string& tag,
	      const std::string& bw,
	      const std::string& fw)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("%s__%s__%s",bw.c_str(),fw.c_str(),tag.c_str());
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  cout<<"List of corr:"<<endl;
	  for(const auto& [key,val] : rawData)
	    cout<<key<<endl;
	  
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	for(size_t t=0;t<T;t++)
	  {
	    union
	    {
	      complex<double> c{};
	      double d[2];
	    };
	    
	    for(size_t ri=0;ri<2;ri++)
	      d[ri]=v[iConf][0][t+T*ri];
	    
	    res[idx({t,iConf})]+=c;
	  }
      
      return res;
    };
  
  auto get=
    [&idx,
     &getRaw,
     &nConfs](const auto& rawData,
	      const std::string& tag,
	      const std::string& bw,
	      const std::string& fw,
	      const bool& ri)
    {
      djvec_t res(T);
      
      const auto d=
	getRaw(rawData,tag,bw,fw);
      
      for(size_t t=0;t<T;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
			 {
			   const double r=d[idx({t,iConf})].real();
			   const double i=d[idx({t,iConf})].imag();
			   res[t][iClust]+=weight*(ri?i:r);
		       });
	}
      res.clusterize(((double)nConfs/njacks));
      
      return res;
    };
  
  const auto getMPCAC=
    [&rawData,
     &get](const std::string& i)
    {
      const djvec_t AP5P5=get(rawData,"P5P5",i,i,0);
      const djvec_t AV0P5=get(rawData,"V0P5",i,i,1);
      
      // effective_mass
      return (forward_derivative(AV0P5.symmetrized(-1))/AP5P5.symmetrized());
    };
  
  effective_mass(get(rawData,"P5P5","H0","H0",0).symmetrized()).ave_err().write("plots/P5P5_run0.xmg");
  effective_mass(get(rawData,"P5P5","H1","H1",0).symmetrized()).ave_err().write("plots/P5P5_run1.xmg");
  
  const djvec_t M0=getMPCAC("H0");
  M0.ave_err().write("plots/mPCAC_run0.xmg");
  const djvec_t M1=getMPCAC("H1");
  M1.ave_err().write("plots/mPCAC_run1.xmg");
  
  // const double k0=0.1394267;
  // const double k1=0.1394;
  // const double k0=0.15;
  // const double k1=0.151;
  // const double k0=0.17;
  // const double k1=0.171;
  // const double k0=0.168331;
  // const double k1=0.168431;
    
  // const double m0=1/k0;
  // const double m1=1/k1;
  // const djvec_t dMPCAC=((M1-M0)/(m1-m0));
  // dMPCAC.ave_err().write("plots/dMPCAC.xmg");
  
  // const djvec_t dm=M0/dMPCAC;
  // const djvec_t k=1/(m0-dm);
  // const djack_t kval=constant_fit(k,20,25,"plots/k.xmg");
  
  // cout<<"k: "<<kval.ave_err()<<endl;
  
  return 0;
}
