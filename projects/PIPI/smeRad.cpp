#include "common.hpp"

int main()
{
  njacks=50;
  
  const std::string dataPath="out";
  const std::vector<std::string> confs=getConfs("confstList.dat",dataPath);
  const size_t nConfs=confs.size();
  
  const auto rawData=
     getRaw("raw.dat",
	    "mes_contr_2Pi",
	    {""},
	    T,
	    dataPath,
	    confs);
  
  auto getRaw=
    [&rawData,
     &nConfs](const std::string& tag)
    {
      const index_t idx({{"t",T},{"conf",nConfs}});
      
      vector<double> res(idx.max());
      
      const auto _v=
	rawData.find(tag);
      if(_v==rawData.end())
	{
	  cout<<"List of corr:"<<endl;
	  for(const auto& [key,val] : rawData)
	    cout<<key<<endl;
	  
	  CRASH("Unable to find %s",tag.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	for(size_t t=0;t<T;t++)
	    res[idx({t,iConf})]=v[iConf][0][t];
      
      return [res,
	      idx](const size_t& t,
		   const size_t& iConf)
      {
	return res[idx({t,iConf})];
      };
    };
  
  // {
  //   const size_t iConf=1;
  //   cout<<confs[iConf]<<endl;
  //   const InterpDef& bSo = interpDef[0];
  //   const InterpDef& bSi = interpDef[0];
    
  //   const std::string& so{bSo.id};
  //   const std::string& si2{bSi.rap[1]};
  //   cout<<so<<" "<<si2<<endl;
    
  //   const auto _A=getRawDirect(so,si2);
  //   for(size_t t=0;t<T;t++)
  //     cout<<_A(0,t,iConf)/(L*L*L)<<endl;
  // }
  
  auto get=
    [&](const string& tag)
    {
      djvec_t res(T);
      
      const auto d=
	getRaw(tag);
      
      for(size_t t=0;t<T;t++)
	jackknivesFill(nConfs,
		       [&](const size_t& iConf,
			   const size_t& iClust,
			   const double& weight)
		       {
			 res[t][iClust]+=weight*d(t,iConf);
		       });
      res.clusterize(((double)nConfs/njacks)).symmetrize();
      
      return res;
    };
  
  for(const auto& [key,val] : rawData)
    {
      const djvec_t k=get(key);
      k.ave_err().write("plots/"+key+".xmg");
      effective_mass(k).ave_err().write("plots/eff_"+key+".xmg");
    }
  
  /////////////////////////////////////////////////////////////////
  
  std::vector<djvec_t> c;
  c.push_back(get("P3_SR0_M3_0__M3_SR0_P3_0__P5P5"));
  c.push_back(get("P3_SMM3_SR0_M3_0__M3_SMP3_SR0_P3_0__P5P5"));
  c.push_back(get("P3_SMM3_SR0_M3_0__M3_SMP3_SR0_P3_0__P5P5"));
  c.push_back(get("P3_SMM3_SR0_SMP3_M3_0__M3_SMP3_SR0_SMM3_P3_0__P5P5"));
  
  const size_t t0=6;
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  tie(eig,recastEigvec,ignore)=gevp(c,t0);
  
  for(size_t i=0;i<2;i++)
    effective_mass(eig[i]).ave_err().write("plots/gevp"+to_string(i)+".xmg");
  
  return 0;
}
