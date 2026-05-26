#include "common.hpp"

void computeDirect()
{
  const std::string directDataPath="outDir";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",directDataPath,"finished");
  const size_t nConfs=confs.size();
  cout<<"NConfs: "<<nConfs<<endl;
  
  const auto rawData=
    getRaw("rawDirect.dat",
	   "mes_contr_dir",
	   {""},
	   T,
	   directDataPath,
	   confs,
	   {"P5P5"});
  
  const size_t nHits=rawData.begin()->second.front().size();
  cout<<"NHits: "<<nHits<<endl;
   
  auto getRawDirect=
    [&rawData,
     &nHits,
     &nConfs](const std::string& bw,
	      const std::string& fw)
    {
      const index_t idx({{"t",T},{"conf",nConfs}});
      
      vector<double> res(idx.max());
      vector<double> par(idx.max());
      vector<double> sin(idx.max());
      
      const string what=
	combine("%s__%s,__P5P5",bw.c_str(),fw.c_str());
      
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
      
      const size_t nHitsToUse=nHits;
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	for(size_t t=0;t<T;t++)
	  {
	    double sumNorm{};
	    double sumProd{};
	    complex<double> sumCompl{};
	    
	    for(size_t iHit=0;iHit<nHitsToUse;iHit++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		sumProd+=(c*c).real();
		sumNorm+=norm(c);
		sumCompl+=c;
	      }
	    
	    res[idx({t,iConf})]=L*L*L*(norm(sumCompl)-sumNorm)/(nHitsToUse*(nHitsToUse-1));
	    par[idx({t,iConf})]=L*L*L*((sumCompl*sumCompl).real()-sumProd)/(nHitsToUse*(nHitsToUse-1));
	    sin[idx({t,iConf})]=sumCompl.real()/nHitsToUse;
	  }
      
      return [data=std::array<std::vector<double>,3>{res,par,sin},
	      idx](const size_t& i,
		   const size_t& t,
		   const size_t& iConf)
      {
	return data[i][idx({t,iConf})];
      };
    };
  
  auto getDirect=
    [&](const string& bw,
	const string& fw)
    {
      std::array<djvec_t,3> res{djvec_t(T),djvec_t(T),djvec_t(T)};
      
      const auto d=
	getRawDirect(bw,fw);
      
      // for(size_t iConf=0;iConf<nConfs;iConf++)
      // 	printf("%s %s %s %.16lg\n",confs[iConf].c_str(),bw.c_str(),fw.c_str(),d(0,1,iConf)/L/L/L);
      // cout<<endl;
      
      for(size_t i=0;i<3;i++)
	{
	  for(size_t t=0;t<T;t++)
	    jackknivesFill(nConfs,
			   [&](const size_t& iConf,
			       const size_t& iClust,
			       const double& weight)
			   {
			     res[i][t][iClust]+=weight*d(i,t,iConf);
			   });
	  res[i].clusterize(((double)nConfs/njacks)).symmetrize();
	}
      
      return res;
    };
  
  enum{DIR,PAR,SIN};
  
  const auto _A=getDirect("bw1","bw1");
  const auto _B=getDirect("bw1","fw1");
  const auto _C=getDirect("bw2","fw2");
  
  grace_file_t dir("plots/dir.xmg");
  dir.set_color_scheme({grace::RED,grace::BLACK,grace::BLUE,grace::GREEN4});
  dir.new_data_set();
  // dir.write_vec_ave_err((2*effective_mass(_A[DIR])).ave_err());
  // dir.set_legend("2Rest");
  dir.write_vec_ave_err(effective_mass(_B[DIR]).ave_err());
  dir.set_legend("Dir");
  dir.write_vec_ave_err(effective_mass(_B[DIR]-_C[DIR]).ave_err());
  dir.set_legend("Dir_proj");
  dir.write_vec_ave_err((2*sqrt(sqr(effective_mass(_A[SIN]))+3*sqr(2*M_PI/L))).ave_err());
  dir.set_legend("2Motion");
  dir.write_vec_ave_err((2*effective_mass(_B[SIN])).ave_err());
  dir.set_legend("2Motion");
}

void computeBox()
{
  const int tMaxBox{25};
  
  int nConfs;
  
  auto getter=
    [&nConfs](const std::string& l)
    {
      const std::string dataPath="outBox";
      const std::vector<std::string> confs=getConfs("confsBoxList"+l+".dat",dataPath,"finished");
      
      nConfs=confs.size();
      
      return
	getRaw("rawBox"+l+".dat",
	       "mes_contr_box"+l,
	       {""},
	       tMaxBox,
	       dataPath,
	       confs,
	       {"S0S0"});
    };
  
  auto rawDataA=getter("A");
  auto rawDataB=getter("B");
  
  auto getRawBox=
    [nHits=1,
     &nConfs](const auto& rawData,
	      const std::string& bw,
	      const std::string& fw)
    {
      const index_t idx({{"t",T},{"conf",nConfs}});
      
      vector<double> res(idx.max());
      
      const string what=
	combine("%s__%s,__S0S0",bw.c_str(),fw.c_str());
      
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
    };
  
  getRawBox(rawDataA,"bw","fwA");
  getRawBox(rawDataB,"bw","fwB");
  
  
  // const size_t nHits=rawData.begin()->second.front().size();
      
  //     const index_t idx({{"hit",nHits},{"tMax",tMaxBox},{"conf",nConfs}});
  
  // auto getRawBox=
  //   [&rawData,
  //    &nHits,
  //    &idx,
  //    &nConfs](const std::string& a,
  // 	      const std::string& b)
  //   {
  //     vector<complex<double>> res(idx.max());
      
  //     const string what=
  // 	combine("%s__%s,__P5P5",a.c_str(),b.c_str());
  //     // cout<<"Searching for "<<what<<endl;
      
  //     const auto _v=
  // 	rawData.find(what);
  //     if(_v==rawData.end())
  // 	CRASH("Unable to find %s",what.c_str());
      
  //     const auto& v=_v->second;
      
  //     for(size_t iHit=0;iHit<nHits;iHit++)
  // 	{
  // 	  for(size_t iConf=0;iConf<nConfs;iConf++)
  // 	    for(size_t t=0;t<tMaxBox;t++)
  // 	      {
  // 		union
  // 		{
  // 		  complex<double> c{};
  // 		  double d[2];
  // 		};
		
  // 		for(size_t ri=0;ri<2;ri++)
  // 		  d[ri]=v[iConf][iHit][t+T*ri];
		
  // 		res[idx({iHit,t,iConf})]+=c;
  // 	      }
  // 	}
      
  //     return res;
  //   };
  
  // auto getBox=
  //   [&]()
  //   {
  //     djvec_t res(tMaxBox);
      
  //     const auto d=
  // 	getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
      
  //     for(size_t t=0;t<tMaxBox;t++)
  // 	{
  // 	  jackknivesFill(nConfs,
  // 			 [&](const size_t& iConf,
  // 			     const size_t& iClust,
  // 			     const double& weight)
  // 		       {
  // 			 double o=0;
  // 			 for(size_t iHit=0;iHit<nHits;iHit++)
  // 			   o+=d[idx({iHit,t,iConf})].real();
  // 			 o/=nHits;
  // 			 res[t][iClust]+=weight*o;
  // 		       });
  // 	}
  //     res.clusterize(((double)nConfs/njacks));
      
  //     return res;
  //   };
  
}

int main()
{
  set_njacks(50);
  
  computeDirect();
  
  computeBox();
  
  return 0;
}
