#include "common.hpp"

djvec_t computeDirect()
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
  dir.set_color_scheme({grace::RED,grace::BLACK,grace::BLUE,grace::GREEN4,grace::ORANGE});
  dir.new_data_set();
  // dir.write_vec_ave_err((2*effective_mass(_A[DIR])).ave_err());
  // dir.set_legend("2Rest");
  dir.write_vec_ave_err(effective_mass(_B[DIR]).ave_err());
  dir.set_legend("Dir");
  dir.write_vec_ave_err(effective_mass(_B[DIR]-_C[DIR]).ave_err());
  dir.set_legend("Dir_proj");
  dir.write_vec_ave_err(effective_mass(_B[DIR]+_C[DIR]).ave_err());
  dir.set_legend("Dir_wr_proj");
  dir.write_vec_ave_err((2*sqrt(sqr(effective_mass(_A[SIN]))+3*sqr(2*M_PI/L))).ave_err());
  dir.set_legend("2Motion");
  dir.write_vec_ave_err((2*effective_mass(_B[SIN])).ave_err());
  dir.set_legend("2Rest");
  dir.write_vec_ave_err((2*effective_mass(_A[SIN])).ave_err());

  grace_file_t inte("plots/dirInte.xmg");
  inte.write_vec_ave_err((_B[DIR]-_C[DIR]).ave_err());
  inte.write_vec_ave_err((_B[DIR]-_C[DIR]-sqr(_B[SIN])*L*L*L).ave_err());
  
  return _B[DIR]-_C[DIR];
}

djvec_t computeBox()
{
  const int tMaxBox{25};
  
  size_t nConfs;
  
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
  
  const index_t idx({{"t",T},{"conf",nConfs}});
  
  auto getRawBox=
    [&idx,
     &nConfs](const auto& rawData,
	      const std::string& bw,
	      const std::string& fw)
    {
      vector<complex<double>> res(idx.max());
      
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
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	for(size_t t=0;t<tMaxBox;t++)
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
  
  auto getBox=
    [&idx,
     &getRawBox,
     &nConfs](const auto& rawData,
	      const std::string& bw,
	      const std::string& fw)
    {
      djvec_t res(tMaxBox);
      
      const auto d=
	getRawBox(rawData,bw,fw);
      
      for(size_t t=0;t<tMaxBox;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
			 {
			   res[t][iClust]+=weight*d[idx({t,iConf})].real();
		       });
	}
      res.clusterize(((double)nConfs/njacks));
      
      return res;
    };
  
  const djvec_t A=getBox(rawDataA,"bw","fwA");
  const djvec_t B=getBox(rawDataB,"bw","fwB");
  
  A.ave_err().write("plots/A.xmg");
  B.ave_err().write("plots/B.xmg");
  
  effective_mass(A).ave_err().write("plots/effA.xmg");
  effective_mass(B).ave_err().write("plots/effB.xmg");
  
  return B-A;
}

int main()
{
  set_njacks(50);
  
  grace_file_t totPlot("plots/tot.xmg");
  
  const djvec_t dir=computeDirect();
  totPlot.write_vec_ave_err(effective_mass(dir).ave_err());
  
  const djvec_t box=computeBox();
  
  const djvec_t tot=dir-2*box;
  totPlot.write_vec_ave_err(effective_mass(tot,T/2,1).ave_err());
  const djvec_t tot2=dir+2*box;
  totPlot.write_vec_ave_err(effective_mass(tot2,T/2,1).ave_err());
  
  return 0;
}
