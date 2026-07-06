#include "common.hpp"

const size_t tMaxBox=26;
const double aGeVInv=0.07948/0.197;;

struct InterpDef
{
  std::string rep;
  
  std::string id;
  
  std::vector<string> rap;
};

const std::vector<InterpDef> interpDef{
  {"001","0",{"PZ","MZ"}},
  {"011","1",{"0P11","0M11"}},
  {"111","2",{"P111","M111"}},
  {"002","3",{"PZ2","MZ2"}},
  {"012","4",{"P012","M012"}}
};

const size_t nOp=interpDef.size();

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeBox(const index_t& idOut)
{
  const std::string boxCorrPath="boxCorr";
  const std::vector<std::string> confs=getConfs("confsBoxList.dat",boxCorrPath,"finished");
  
  const size_t nConfs=confs.size();
  
  const auto rawData=
    getRaw("rawBox.dat",
	   "mes_contr_box",
	   {""},
	   tMaxBox,
	   boxCorrPath,
	   confs);
  const size_t nHits=rawData.begin()->second.front().size();
  
  const index_t idx({{"hit",nHits},{"tMax",tMaxBox},{"conf",nConfs}});
  
  auto getRawBox=
    [&rawData,
     &nHits,
     &idx,
     &nConfs](const std::string& a,
	      const std::string& b)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("%s__%s,__P5P5",a.c_str(),b.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  for(const auto& [tag,vale] : rawData)
	    cout<<" "<<tag<<endl;
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<tMaxBox;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		res[idx({iHit,t,iConf})]+=c;
	      }
	}
      
      return res;
    };
  
  auto getBox=
    [&](const string& mso1,
	const string& mso2,
	const string& msi1,
	const string& msi2)
    {
      djvec_t res(tMaxBox);
      
      const auto d=
	// getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
	getRawBox("bw"+mso1,"fw"+msi2);
      
      for(size_t t=0;t<tMaxBox;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
		       {
			 double o=0;
			 for(size_t iHit=0;iHit<nHits;iHit++)
			   o+=d[idx({iHit,t,iConf})].real();
			 o/=nHits;
			 res[t][iClust]+=weight*o;
		       });
	}
      res.clusterize(((double)nConfs/njacks));
      
      return res;
    };
  
  std::vector<djvec_t> res(idOut.max());
  
  for(size_t iBSo=0;iBSo<interpDef.size();iBSo++)
    for(size_t iBSi=0;iBSi<interpDef.size();iBSi++)
      {
	const InterpDef& bSo{interpDef[iBSo]};
	const InterpDef& bSi{interpDef[iBSi]};
	
	const std::string& repSo{bSo.rep};
	const std::string& repSi{bSi.rep};
	const std::string& so1{bSo.rap[0]};
	const std::string& so2{bSo.rap[1]};
	const std::string& si{bSi.id};
	const djvec_t A=getBox(so1,so2,"M"+si,"P"+si);
	const djvec_t B=getBox(so2,so1,"M"+si,"P"+si);
	
	// effective_mass(g,T/2).ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
	// effective_mass(h,T/2).ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	A.ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
	B.ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	const djvec_t C=A-B;
	C.ave_err().write("plots/C_"+repSi+"_"+repSo+".xmg");
	
	res[idOut({iBSo,iBSi})]=C;
      }
  
  return res;
}

auto box()
{
  index_t idOut({{"bSo",interpDef.size()},{"bSi",interpDef.size()}});
  
  return [data=computeOrLoad(idOut,"box.dat",computeBox),
	  idOut](const size_t& iBso,
		 const size_t& iBsi) -> const djvec_t&
  {
    return data[idOut({iBso,iBsi})];
  };
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeCurrent(const index_t&)
{
  const std::string corrPath="directCorr";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",corrPath,"finished");
  
  const size_t nConfs=confs.size();
  
  const auto rawData=
    getRaw("rawCur.dat",
	   "mes_contr_current",
	   {""},
	   T,
	   corrPath,
	   confs,
	   {"V1V1","V2V2","V3V3"});
  const size_t nHits=rawData.begin()->second.front().size();
  
  const index_t idx({{"hit",nHits},{"tMax",T},{"conf",nConfs}});
  
  auto getRawCur=
    [&rawData,
     &nHits,
     &idx,
     &nConfs](const std::string& tag)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("propR1__propR0,__%s",tag.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  for(const auto& [tag,vale] : rawData)
	    cout<<" "<<tag<<endl;
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<T;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		res[idx({iHit,t,iConf})]+=c;
	      }
	}
      
      return res;
    };
  
  const auto d=
    // getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
    (getRawCur("V1V1")+getRawCur("V2V2")+getRawCur("V3V3"))*(1.0/3);
  
  auto comb=
    [&](const auto& f)
    {
      djvec_t res(T);
      for(size_t t=0;t<T;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
			 {
			   double o=0;
			   for(size_t iHit=0;iHit<nHits;iHit++)
			     o+=d[idx({iHit,t,iConf})].real();
			   o/=nHits;
			   res[t][iClust]+=weight*o;
			 });
	}
      res.clusterize(((double)nConfs/njacks)).symmetrize();
      
      return res;
    };
    
  return {comb(d)};
}

auto current()
{
  index_t idOut({{"dum",1}});
  
  return [data=computeOrLoad(idOut,"cur.dat",computeCurrent)](const size_t i) -> const djvec_t&
  {
    return data[i];
  };
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computePion(const index_t&)
{
  const std::string corrPath="directCorr";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",corrPath,"finished");
  
  const size_t nConfs=confs.size();
  
  const auto rawData=
    getRaw("rawPi.dat",
	   "mes_contr_pi",
	   {""},
	   T,
	   corrPath,
	   confs,
	   {"P5P5"});
  const size_t nHits=rawData.begin()->second.front().size();
  
  const index_t idx({{"hit",nHits},{"tMax",T},{"conf",nConfs}});
  
  auto getRawPi=
    [&rawData,
     &nHits,
     &idx,
     &nConfs](const std::string& tag)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("propR0__propR0,__%s",tag.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  for(const auto& [tag,vale] : rawData)
	    cout<<" "<<tag<<endl;
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<T;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		res[idx({iHit,t,iConf})]+=c;
	      }
	}
      
      return res;
    };
  
  const auto d=
    getRawPi("P5P5");
  
  auto comb=
    [&](const auto& f)
    {
      djvec_t res(T);
      for(size_t t=0;t<T;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
			 {
			   double o=0;
			   for(size_t iHit=0;iHit<nHits;iHit++)
			     o+=d[idx({iHit,t,iConf})].real();
			   o/=nHits;
			   res[t][iClust]+=weight*o;
			 });
	}
      res.clusterize(((double)nConfs/njacks)).symmetrize();
      
      return res;
    };
    
  return {comb(d)};
}

auto pion()
{
  index_t idOut({{"dum",1}});
  
  return [data=computeOrLoad(idOut,"pi.dat",computePion)](const size_t i) -> const djvec_t&
  {
    return data[i];
  };
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeTri(const index_t& idOut)
{
  const std::string corrPath="triCorr";
  const std::vector<std::string> confs=getConfs("confsTriList.dat",corrPath,"finished");
  
  const size_t nConfs=confs.size();
  
  map<string,vector<vector<vector<double>>>> rawData;
  for(const auto& [mesFile,rawFile,what] : std::vector<std::array<string,3>>{
      {"mes_contr_tri","rawTri","V3P5"},
      {"mes_contr_triT","rawTriT","T3P5"},
      {"mes_contr_triSme","rawTriSme","V3P5"},
      {"mes_contr_triTSme","rawTriTSme","T3P5"}})
    for(const auto& it : getRaw(rawFile,
				mesFile,
				{""},
				T,
				corrPath,
				confs,
				{what}))
      rawData.insert(it);
  
  const size_t nHits=rawData.begin()->second.front().size();
  
  const index_t idx({{"hit",nHits},{"tMax",T},{"conf",nConfs}});
  
  auto getRawTri=
    [&rawData,
     &nHits,
     &idx,
     &nConfs](const std::string& p,
	      const std::string& maybeSm,
	      const std::string& tag)
    {
      vector<complex<double>> res(idx.max());
      
      const string what=
	combine("bw%s__sm%s_prop,__%s",p.c_str(),maybeSm.c_str(),tag.c_str());
      // cout<<"Searching for "<<what<<endl;
      
      const auto _v=
	rawData.find(what);
      if(_v==rawData.end())
	{
	  for(const auto& [tag,vale] : rawData)
	    cout<<" "<<tag<<endl;
	  CRASH("Unable to find %s",what.c_str());
	}
      
      const auto& v=_v->second;
      
      for(size_t iHit=0;iHit<nHits;iHit++)
	{
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    for(size_t t=0;t<T;t++)
	      {
		union
		{
		  complex<double> c{};
		  double d[2];
		};
		
		for(size_t ri=0;ri<2;ri++)
		  d[ri]=v[iConf][iHit][t+T*ri];
		
		res[idx({iHit,t,iConf})]+=c;
	      }
	}
      
      return res;
    };
  
  auto getTri=
    [&](const string& repSo,
	const string& p,
	const string& what,
	const string& maybeSme,
	const int& parity)
    {
      djvec_t res(T);
      
      vector<complex<double>> d =
          // getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
	getRawTri(p,maybeSme,what);
      
      for(size_t t=0;t<T;t++)
	{
	  jackknivesFill(nConfs,
			 [&](const size_t& iConf,
			     const size_t& iClust,
			     const double& weight)
		       {
			 double o=0;
			 // const int mu=2;
			 for(size_t iHit=0;iHit<nHits;iHit++)
			   // for(int mu=0;mu<3;mu++)
			     // if(repSo[mu]!='0')
			       o+=d[idx({iHit,t,iConf})].imag();
			 res[t][iClust]+=weight*o/nHits;
		       });
	}
      res.clusterize(((double)nConfs/njacks)).symmetrize(parity);
      
      return res;
    };
  
  std::vector<djvec_t> res(idOut.max());
  
  for(size_t iBSo=0;iBSo<interpDef.size();iBSo++)
    for(size_t iSm=0;iSm<2;iSm++)
      for(size_t iVT=0;iVT<2;iVT++)
	{
	  const InterpDef& bSo{interpDef[iBSo]};
	  const std::string& repSo{bSo.rep};
	  const std::string& so1{bSo.rap[0]};
	  const std::string& so2{bSo.rap[1]};
	  
	  auto g=
	    [&repSo,
	     &iSm,
	     &iVT,
	     &getTri](const string& so)
	    {
	      return getTri(repSo,so,(iVT==0)?"V3P5":"T3P5",(iSm==0)?"":"_sm",(iVT==0)?1:-1);
	    };
	  const djvec_t A=g(so1);
	  const djvec_t B=g(so2);
	  
	  // effective_mass(g,T/2).ave_err().write("plots/A_"+repSi+"_"+repSo+".xmg");
	  // effective_mass(h,T/2).ave_err().write("plots/B_"+repSi+"_"+repSo+".xmg");
	  A.ave_err().write("plots/triA_"+repSo+".xmg");
	  B.ave_err().write("plots/triB_"+repSo+".xmg");
	  const djvec_t C=A-B;
	  C.ave_err().write("plots/triC_"+repSo+".xmg");
	  
	  res[idOut({iBSo,iSm,iVT})]=C/sqrt(2);
	  // getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
	}
  
  return res;
}

auto triangle()
{
  index_t idOut{{{{"bSo",interpDef.size()},{"sm",2},{"vt",2}}}};
  
  return [data=computeOrLoad(idOut,"triangle.dat",computeTri),
	  idOut](const size_t& iOp,
		 const size_t& iSm,
		 const size_t& iVt) -> const djvec_t&
  {
    return data[idOut({iOp,iSm,iVt})];
  };
}

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeDirect(const index_t& idOut)
{
  const std::string directDataPath="directCorr";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",directDataPath,"finished");
  const size_t nConfs=confs.size();
  
  const auto rawData=
     getRaw("rawDirect.dat",
	    "mes_contr_dir",
	    {""},
	    T,
	    directDataPath,
	    confs);
  
  const size_t nHits=rawData.begin()->second.front().size();
  cout<<"NHits: "<<nHits<<endl;
  
  auto getRawDirect=
    [&rawData,
     &nHits,
     &nConfs](const std::string& mso,
	      const std::string& msi)
    {
      const index_t idx({{"t",T},{"conf",nConfs}});
      
      vector<double> res(idx.max());
      vector<double> par(idx.max());
      vector<double> sin(idx.max());
      
      const string what=
	combine("bw%s_%s__fw%s_%s,__P5P5",msi.c_str(),mso.c_str(),msi.c_str(),mso.c_str());
      
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
  
  auto getDirect=
    [&](const string& mso,
	const string& msi)
    {
      std::array<djvec_t,3> res{djvec_t(T),djvec_t(T),djvec_t(T)};
      
      const auto d=
	getRawDirect(mso,msi);
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	printf("%s %s %s %.16lg\n",confs[iConf].c_str(),mso.c_str(),msi.c_str(),d(0,1,iConf)/L/L/L);
      cout<<endl;
      
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
  
  std::vector<djvec_t> res(idOut.max());
  
  for(size_t iBSo=0;iBSo<interpDef.size();iBSo++)
    for(size_t iBSi=0;iBSi<interpDef.size();iBSi++)
      {
	const InterpDef& bSo{interpDef[iBSo]};
	const InterpDef& bSi{interpDef[iBSi]};
	
	const std::string& so{bSo.rap[0]};
	const std::string& si1{bSi.rap[0]};
	const std::string& si2{bSi.rap[1]};
	
	const auto _A=getDirect(so,si1);
	const auto _B=getDirect(so,si2);
	
	for(size_t iCombo=0;iCombo<3;iCombo++)
	  res[idOut({iBSo,iBSi,iCombo})]=_A[iCombo]-_B[iCombo];
	
	const djvec_t A=_A[DIR];
	const djvec_t B=_B[DIR];
	const djvec_t D=_A[PAR];
	const djvec_t E=_B[PAR];
	const djvec_t G=_A[SIN];
	const djvec_t H=_B[SIN];
	
	const std::string& repSo{bSo.rep};
	const std::string& repSi{bSi.rep};
	A.ave_err().write("plots/dA_"+repSi+"_"+repSo+".xmg");
	B.ave_err().write("plots/dB_"+repSi+"_"+repSo+".xmg");
	const djvec_t C=A-B;
	C.ave_err().write("plots/dC_"+repSi+"_"+repSo+".xmg");
	
	(effective_mass(_B[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/effDinte_"+repSi+"_"+repSo+".xmg");
	(effective_mass(_B[DIR]-_A[DIR])-2*effective_mass(_B[SIN])).ave_err().write("plots/effDinte2_"+repSi+"_"+repSo+".xmg");
	
	effective_mass(A).ave_err().write("plots/effDA_"+repSi+"_"+repSo+".xmg");
	effective_mass(B).ave_err().write("plots/effDB_"+repSi+"_"+repSo+".xmg");
	effective_mass(D).ave_err().write("plots/effDD_"+repSi+"_"+repSo+".xmg");
	effective_mass(E).ave_err().write("plots/effDE_"+repSi+"_"+repSo+".xmg");
	(2*effective_mass(G)).ave_err().write("plots/effDG_"+repSi+"_"+repSo+".xmg");
	(2*effective_mass(H)).ave_err().write("plots/effDH_"+repSi+"_"+repSo+".xmg");
	effective_mass(A-B).ave_err().write("plots/effDC_"+repSi+"_"+repSo+".xmg");
	effective_mass(H-G).ave_err().write("plots/effDI_"+repSi+"_"+repSo+".xmg");
    }
  
  return res;
}

auto direct()
{
  index_t idOut{{{"bSo",interpDef.size()},{"bSi",interpDef.size()},{"combo",3}}};
  
  return [data=computeOrLoad(idOut,"direct.dat",computeDirect),
	  idOut](const size_t& iBso,
		 const size_t& iBsi,
		 const size_t& iCombo) -> const djvec_t&
  {
    return data[idOut({iBso,iBsi,iCombo})];
  };
}

int main()
{
  njacks=50;
  
  const auto b=box();
  
  const auto d=direct();
  
  const auto jj=current();
  
  const auto pi=pion();
  
  const auto tri=triangle();
  
  const size_t nOpToUse=nOp;
  
  index_t iC({{"so",nOpToUse+1},{"si",nOpToUse+1}});
  std::vector<djvec_t> c((nOpToUse+1)*(nOpToUse+1));
  
  for(size_t ibSo=0;ibSo<nOpToUse;ibSo++)
    {
      (tri(ibSo,1/*ism*/,/*iVt*/ 0)/tri(ibSo,0/*ism*/,/*iVt*/ 0)).ave_err().write(combine("plots/rat_SM_noSM%zu.xmg",ibSo));
      (-tri(ibSo,0/*ism*/,/*iVt*/ 1)/tri(ibSo,0/*ism*/,/*iVt*/ 0)).ave_err().write(combine("plots/rat_T_V_%zu.xmg",ibSo));
    }
  
  for(size_t ibSo=0;ibSo<nOpToUse;ibSo++)
    {
      for(size_t ibSi=0;ibSi<nOpToUse;ibSi++)
	{
	  const djvec_t& box=b(ibSo,ibSi);
	  const djvec_t& dir=-d(ibSo,ibSi,0); //so it is
	  box.ave_err().write(combine("plots/box_%zu_%zu.xmg",ibSo,ibSi));
	  dir.ave_err().write(combine("plots/dir_%zu_%zu.xmg",ibSo,ibSi));
	  const djvec_t cmb=dir-2*box;
	  
	  c[iC({ibSo+1,ibSi+1})]=cmb;
	  cmb.ave_err().write(combine("plots/cmb_%zu_%zu.xmg",ibSo,ibSi));
	}
      
      c[iC({ibSo+1,0})]=
	c[iC({0,ibSo+1})]=
	tri(ibSo,0,0);
    }
  
  c[0]=jj(0);
  
  for(size_t i=0;i<nOpToUse+1;i++)
    for(size_t j=0;j<nOpToUse+1;j++)
      {
	djvec_t& a=
	  c[iC({i,j})];
	djvec_t& b=
	  c[iC({j,i})];
	
	a=b=(a+b).subset(0,tMaxBox-1)/2;
      }
  
  effective_mass(jj(0)).ave_err().write("plots/jj.xmg");
  
  const djack_t mPi=constant_fit(effective_mass(pi(0)),15,20,"plots/pi.xmg");
  
  const size_t t0=7;
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  for(size_t i=0;i<nOpToUse+1;i++)
    {
      for(size_t j=0;j<nOpToUse+1;j++)
	cout<<c[iC({i,j})][10].ave_err()<<"     ";
      
      cout<<endl;
    }
  
  cout<<"VV: "<<combine("%.16lg\n",c[0][10].ave())<<endl;
  cout<<"tri: "<<combine("%.16lg\n",c[1][10].ave())<<endl;
  
  tie(eig,recastEigvec,ignore)=gevp(c,t0);
  
  // typedef Matrix<double,Dynamic,Dynamic> Matr;
  // GeneralizedEigenSolver<Matr> ges;
  // Matr aM(nOpToUse+1,nOpToUse+1);
  // Matr bM(nOpToUse+1,nOpToUse+1);
  // vector<djvec_t> eig2(nOpToUse+1,djvec_t(tMaxBox));
  // for(size_t ijack=0;ijack<=njacks;ijack++)
  //   {
  //     //Fill the rhs matrix
  //     for(size_t i=0;i<nOpToUse+1;i++)
  // 	for(size_t j=0;j<nOpToUse+1;j++)
  // 	  bM(i,j)=c[iC({i,j})][t0][ijack];
      
  //     for(size_t t=0;t<tMaxBox;t++)
  // 	{
  // 	  for(size_t i=0;i<nOpToUse+1;i++)
  // 	    for(size_t j=0;j<nOpToUse+1;j++)
  // 	      aM(i,j)=c[iC({i,j})][t][ijack];
	  
	  
  // 	  ges.compute(aM,bM,1);
	  
  // 	  for(size_t i=0;i<nOpToUse+1;i++)
  // 	    eig2[i][t][ijack]=ges.eigenvalues()[i].real();
  // 	}
  //   }
  
  vector<double> expSh{0.0019364276101050126,0.016167593494143095,0.028595040062616484,0.04024620996370226,0.04075042586371391,0.04075042586371391};
  
  // auto resort=
  //   [nOpToUse](std::vector<djvec_t>& eig)
  //   {
  for(size_t t=0;t<tMaxBox;t++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      {
	std::vector<double> temp(nOpToUse+1);
	for(size_t iop=0;iop<nOpToUse+1;iop++)
	  temp[iop]=eig[iop][t][ijack];
	std::sort(temp.begin(),temp.end());
	for(size_t iop=0;iop<nOpToUse+1;iop++)
	  eig[iop][t][ijack]=temp[nOpToUse-iop];
      }
// };
  
//   resort(eig);
//   resort(eig2);
  
//   for(size_t ijack=0;ijack<=njacks;ijack++)
//     for(size_t t=0;t<tMaxBox;t++)
//       for(size_t i=0;i<nOpToUse+1;i++)
// 	cout<<ijack<<" jack, "<<t<<" t, "<<i<<" op, "<<eig[i][t][ijack]<<" "<<eig2[i][t][ijack]<<endl;
  
  // for(size_t t=0;t<tMaxBox;t++)
  //   {
  //     for(size_t iop=0;iop<nOpToUse+1;iop++)
  // 	cout<<eig[iop][t].ave()<<" ";
  //     cout<<endl;
  //   }
  
  const std::vector<double> n2={1,2,3,4,5,6};
  
  grace_file_t gevp("plots/gevp.xmg");
  gevp.set_color_scheme({grace::BLACK,grace::BLUE,grace::RED,grace::ORANGE,grace::GREEN4,grace::VIOLET});
  for(size_t iop=0;iop<nOpToUse+1;iop++)
    {
      const djvec_t eigP=effective_mass(eig[iop],T/2);
      gevp.write_vec_ave_err(eigP.ave_err(),gevp.color_scheme[iop%gevp.color_scheme.size()],gevp.symbol_scheme[iop&gevp.symbol_scheme.size()]);
      gevp.set_no_line();
      
      cout<<iop<<" "<<eigP[10].ave_err()<<endl;
      
      const djack_t E=2*sqrt(sqr(mPi)+sqr(2*M_PI/64)*n2[iop])-expSh[iop];
      gevp.write_constant_band(0,tMaxBox,E,gevp.color_scheme[iop]);
    }
  
  grace_file_t fullCompa("plots/fullCompa.xmg");
  fullCompa.set_color_scheme({grace::BLACK,grace::BLUE,grace::RED,grace::ORANGE,grace::GREEN4});
  for(size_t iop=0;iop<nOpToUse;iop++)
    {
      grace_file_t cmp("plots/cmp"+to_string(iop)+".xmg");
      cmp.set_title(interpDef[iop].rep);
      cmp.set_color_scheme({grace::BLACK,grace::BLUE,grace::RED,grace::ORANGE});
      cmp.set_no_line();
      
      const djvec_t nonInter=(2*effective_mass(d(iop,iop,2),T/2));
      cmp.write_vec_ave_err(nonInter.ave_err());
      cmp.set_legend("Non interacting pion");
      cmp.new_data_set();
      cmp.set_no_line();
      
      const djvec_t eigP=effective_mass(eig[iop],T/2);
      cmp.write_vec_ave_err(eigP.ave_err());
      cmp.set_no_line();
      cmp.set_legend("gevp");
      
      const djvec_t full=effective_mass(c[iC({iop+1,iop+1})],T/2);
      cmp.write_vec_ave_err(full.ave_err());
      cmp.set_no_line();
      cmp.set_legend("2Box-direct");
      
      fullCompa.write_vec_ave_err(full.ave_err());
      fullCompa.set_legend(interpDef[iop].rep);
      fullCompa.set_no_line();
      
      const djack_t E=2*sqrt(sqr(mPi)+sqr(2*M_PI/64)*n2[iop])-expSh[iop];
      // const djvec_t expected=(2*effective_mass(d(0,0,2),T/2)-expSh[iop]);
      cmp.write_constant_band(0,tMaxBox,E);
      cmp.set_legend("GS expectation");
      
      grace_file_t shiftPlot("plots/shift"+to_string(iop)+".xmg");
      shiftPlot.write_vec_ave_err((full-nonInter).ave_err());
      shiftPlot.write_line([e=-expSh[iop]](const double&){return e;},0,tMaxBox);
      
      //cout<<expSh[iop]/aGeVInv<<endl;
      
      // /////////////////////////////////////////////////////////////////
      
      grace_file_t gevpCompa("plots/gevp"+to_string(iop)+".xmg");
      gevpCompa.set_title(interpDef[iop].rep);
      gevpCompa.set_color_scheme({grace::BLACK,grace::BLUE,grace::RED,grace::ORANGE});
      gevpCompa.set_no_line();
      
      gevpCompa.write_vec_ave_err(full.ave_err());
      gevpCompa.set_legend("2Box-direct");
      gevpCompa.new_data_set();
      gevpCompa.set_no_line();
      
      const djvec_t eigEffMass=effective_mass(eig[iop],T/2);
      gevpCompa.write_vec_ave_err(eigEffMass.ave_err());
      gevpCompa.set_no_line();
      gevpCompa.set_legend("Eig");
    }
  
  return 0;
}
