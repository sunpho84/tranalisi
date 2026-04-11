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
  c.push_back(get("P3_SMP3_SR0_M3_0__M3_SMM3_SR0_P3_0__P5P5"));
  c.push_back(get("P3_SMP3_SR0_M3_0__M3_SMM3_SR0_P3_0__P5P5"));
  c.push_back(get("P3_SMP3_SR0_SMP3_M3_0__M3_SMM3_SR0_SMM3_P3_0__P5P5"));
  
  const double x=0.3;
  const djvec_t one=c[3]+c[0]*x*x-x*(c[1]+c[2]);
  effective_mass(one).ave_err().write("plots/one.xmg");
  const double y=1.05;
  const djvec_t two=c[0]+c[3]*y*y-y*(c[1]+c[2]);
  effective_mass(two).ave_err().write("plots/two.xmg");
  
  djvec_t ZL(2);
  djvec_t ZS(2);
  djvec_t M(2);
  jack_fit_t jf;
  jf.add_fit_par(ZL[0],"ZL0",0.19,0.01);
  jf.add_fit_par(ZL[1],"ZL1",0.019,0.001);
  jf.add_fit_par(ZS[0],"ZS0",0.23,0.01);
  jf.add_fit_par(ZS[1],"ZS1",0.004,0.001);
  jf.add_fit_par(M[0],"M0",0.18,0.01);
  jf.add_fit_par(M[1],"M1",0.57,0.01);
  
  for(size_t t=6;t<20;t++)
    {
      auto add=
	[&jf,
	 t](const djvec_t& c,
	    const size_t& LS1,
	    const size_t& LS2)
      {
	jf.add_point(c[t],
		     [t,LS1,LS2](const vector<double> &p,size_t iel){return two_pts_corr_fun(p[0+2*LS1]*p[0+2*LS2],p[4],T/2.0,t,1)+two_pts_corr_fun(p[1+2*LS1]*p[1+2*LS2],p[5],T/2.0,t,1);});
      };
      
      add(c[0],0,0);
      add(c[1],1,0);
      add(c[2],1,1);
    }
  
  jf.fit();
  cout<<M.ave_err()<<endl;
  
  
  const size_t t0=6;
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  tie(eig,recastEigvec,ignore)=gevp(c,t0);
  
  for(size_t i=0;i<2;i++)
    effective_mass(eig[i]).ave_err().write("plots/gevp"+to_string(i)+".xmg");
  
  return 0;
}
