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
  {"001","1",{"Pz","Mz"}},
  {"011","2",{"0P11","0M11"}},
  {"111","3",{"P111","M111"}},
  {"002","4",{"Pz2","Mz2"}},
  {"012","5",{"0P12","0M12"}}
};

const size_t nOp=interpDef.size();

/////////////////////////////////////////////////////////////////

std::vector<djvec_t> computeBox(const index_t& idOut)
{
  const std::string boxCorrPath="boxCorr";
  const std::vector<std::string> confs=getConfs("confsBoxList.dat",boxCorrPath);
  
  const size_t nConfs=confs.size();
  
  const auto rawData=
    getRaw("rawBox.dat",
	   "mes_contr_box_src_snk%s_T25",
	   {"1","2","3","4","5"},
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
	CRASH("Unable to find %s",what.c_str());
      
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
	getRawBox("Sr1_"+mso1+"_D0_G5_Sr0_"+mso2+"_0",msi1+"_TH25_Sr1_G5_"+msi2+"_Sr0_0");
      
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

std::vector<djvec_t> computeDirect(const index_t& idOut)
{
  const std::string directDataPath="directCorr";
  const std::vector<std::string> confs=getConfs("confsDirectList.dat",directDataPath);
  const size_t nConfs=confs.size();
  
  const auto rawData=
     getRaw("rawDirect.dat",
	    "mes_contr_direct_r0_P%s",
	    {"0","1","2","3","4","5"},
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
	combine("M0_Sr0_P0_0__%s_Sr0_P%s_0,__P5P5",msi.c_str(),mso.c_str());
      
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
	
	const std::string& so{bSo.id};
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
	(A-B).ave_err().write("plots/dC_"+repSi+"_"+repSo+".xmg");
	
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
  
  const size_t nOpToUse=nOp;
  std::vector<djvec_t> c(nOpToUse*nOpToUse);
  
  for(size_t ibSo=0;ibSo<nOpToUse;ibSo++)
    for(size_t ibSi=0;ibSi<nOpToUse;ibSi++)
      {
	const djvec_t& box=b(ibSo,ibSi);
	const djvec_t& dir=d(ibSo,ibSi,0);
	box.ave_err().write(combine("plots/box_%zu_%zu.xmg",ibSo,ibSi));
	(-dir).ave_err().write(combine("plots/dir_%zu_%zu.xmg",ibSo,ibSi));
	const djvec_t cmb=2*box-dir;
	c[ibSo+nOpToUse*ibSi]=cmb;
	cmb.ave_err().write(combine("plots/cmb_%zu_%zu.xmg",ibSo,ibSi));
	
	if(ibSi<ibSo)
	  c[ibSo+nOpToUse*ibSi]=c[ibSi+nOpToUse*ibSo];
      }
  
  const size_t t0=10;
  
  vector<djvec_t> eig;
  vector<djvec_t> recastEigvec;
  vector<djvec_t> origEigvec;
  
  tie(eig,recastEigvec,ignore)=gevp(c,t0);
  vector<double> expSh{0.0019364276101050126,0.016167593494143095,0.028595040062616484,0.04024620996370226,0.04075042586371391};
  
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
      
      const djvec_t full=effective_mass(c[iop+nOpToUse*iop],T/2);
      cmp.write_vec_ave_err(full.ave_err());
      cmp.set_no_line();
      cmp.set_legend("2Box-direct");
      
      fullCompa.write_vec_ave_err(full.ave_err());
      fullCompa.set_legend(interpDef[iop].rep);
      fullCompa.set_no_line();
      
      const djvec_t expected=(2*effective_mass(d(iop,iop,2),T/2)-expSh[iop]);
      cmp.write_vec_ave_err(expected.ave_err());
      cmp.set_no_line();
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
