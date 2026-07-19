#include "../PIPI/common.hpp"

template <typename A>
djvec_t get(const A& rawData,
	    const std::vector<string>& confs,
	    const std::string& path,
	    const std::string& a,
	    const std::string& b,
	    const std::string& c)
{
  const string what=
    a+"__"+b+",__"+c;
      
  const auto _v=
    rawData.find(what);
  if(_v==rawData.end())
    {
      for(const auto& [tag,vale] : rawData)
	cout<<" "<<tag<<endl;
      CRASH("Unable to find %s",what.c_str());
    }
  
  const auto& v=_v->second;
  
  djvec_t res(T);
  
  jackknivesFill(confs.size(),
		 [&res,
		  &v](const size_t& iConf,
		     const size_t& iClust,
		     const double& weight)
		 {
		   for(size_t t=0;t<T;t++)
		     res[t][iClust]+=weight*v[iConf][0][t];
		 });
  
  res.clusterize(((double)confs.size()/njacks));
  
  return res;
}

auto dum=
  []()
  {
    njacks=50;
    return 1;
  }();

struct EnsPars
{
  int L;
  
  djack_t a;
  
  EnsPars(const int& L,
	  const double& aAve,
	  const double& eErr,
	  const int seed) :
    L(L)
  {
    cout<<smart_print({aAve,eErr})<<endl;
    a.fill_gauss(aAve/0.19731,eErr/0.19731,seed);
  }
};

const std::map<std::string,EnsPars> ensParsList{
  {"A48", {48,0.090700,0.000500,223525}},
  {"B64", {64,0.079480,0.000110,223526}},
  {"C80", {80,0.068190,0.000140,223527}},
  {"D96", {96,0.056850,0.000090,223528}},
  {"E112",{112,0.048920,0.000110,223529}},
};

int main()
{
  const std::string ensName=
    std::filesystem::path{pwd()}.filename();
  
  const auto it=
    ensParsList.find(ensName);
  if(it==ensParsList.end())
    {
      cerr<<"available ensembles:"<<endl;
      for(const auto& [name,_] : ensParsList)
	cout<<" \""<<name<<"\""<<endl;
      CRASH("Unable to find ens %s",ensName.c_str());
    }
  
  const EnsPars& ensPars=
    it->second;

  L=ensPars.L;
  T=L*2;
  const djack_t& a=ensPars.a;
  
  //set_njacks(50);
  
  const std::string path="out";
  const std::vector<std::string> confs=getConfs("",path,"finished");
  
  const size_t nConfs=confs.size();
  cout<<"NConfs: "<<nConfs<<endl;
  
  const auto rawData=
    getRaw("",
	   "mes_contr_c",
	   {""},
	   T,
	   path,
	   confs,
	   {"P5P5","V1V1","V2V2","V3V3","V0P5"});
  
  std::vector<djvec_t> P5P5(7);
  std::vector<djvec_t> VKVK(7,djvec_t(T/2+1));
  for(size_t i=0;i<7;i++)
    {
      P5P5[i]=get(rawData,confs,path,"C1",combine("C%zu",i+1),"P5P5").symmetrized();
      for(int mu=1;mu<=3;mu++)
	VKVK[i]+=get(rawData,confs,path,"C1",combine("C%zu",i+1),combine("V%dV%d",mu,mu)).symmetrized();
      VKVK[i]/=3;
    }
  const double p0=M_PI*1e-3/L;
  // const double p1=M_PI*1e-5/L;
  // const double p2=2*M_PI/L;
  
  const djvec_t m0=effective_mass(P5P5[0]);
  const djvec_t v0=effective_mass(VKVK[0]);
  const djvec_t e1=effective_mass(P5P5[1]+P5P5[2]);
  // const djvec_t e2=effective_mass(c[3]+c[4]);
  // const djvec_t e3=effective_mass(c[5]+c[6]);
  
  const djvec_t c1=(sqr(e1)-sqr(m0))/sqr(p0);
  // const djvec_t c2=(sqr(e2)-sqr(m0))/sqr(p1);
  // const djvec_t c3=(sqr(e3)-sqr(m0))/sqr(p2);
  
  // const djack_t M=constant_fit(m0/ensPars.a,15,18,"plots/m0.xmg");
  // cout<<M.ave_err()<<endl;
  const djack_t V=constant_fit(v0/a,30,35,"plots/v0.xmg");
  cout<<sqr(a).ave()<<" "<<V.ave_err()<<endl;
  c1.ave_err().write("plots/c1.xmg");
  // c2.ave_err().write("plots/c2.xmg");
  // c3.ave_err().write("plots/c3.xmg");
  // const djvec_t A=read_djvec("mes_contr_P5P5_cp",T).symmetrized();
  // const djvec_t B=read_djvec("mes_contr_P5P5_cm",T).symmetrized();
  // const djvec_t C=read_djvec("mes_contr_P5P5_c0",T).symmetrized();
  
  // effective_mass(A).ave_err().write("/tmp/A.xmg");
  // effective_mass(B).ave_err().write("/tmp/B.xmg");
  // effective_mass(C).ave_err().write("/tmp/C.xmg");
  
  // auto e2=
  //   [](const djvec_t& c)
  // {
  //   return sqr(effective_mass(c));
  // };

  // const double p2=sqr(1e-3*M_PI/64);
  // const djvec_t D=(e2(A)+e2(B)-2*e2(C))/(2*p2);
  // D.ave_err().write("/tmp/D.xmg");
  
  // const djvec_t E=D/C;
  // E.ave_err().write("/tmp/E.xmg");
  
  // return 0;
  
  // const std::string dataPath="out";
  // const std::vector<std::string> confs=getConfs("confsDataList.dat",dataPath,"");
  // const size_t nConfs=confs.size();
  // cout<<"NConfs: "<<nConfs<<endl;
  
  // const auto rawData=
  //   getRaw("RUN0.dat",
  // 	   "mes_contr_HH",
  // 	   {""},
  // 	   T,
  // 	   dataPath,
  // 	   confs,
  // 	   {"P5P5","V0P5","V1V1"});
  
  // const index_t idx({{"t",T},{"conf",nConfs}});
  
  // auto getRaw=
  //   [&idx,
  //    &nConfs](const auto& rawData,
  // 	      const std::string& tag,
  // 	      const std::string& bw,
  // 	      const std::string& fw)
  //   {
  //     vector<complex<double>> res(idx.max());
      
  //     const string what=
  // 	combine("%s__%s__%s",bw.c_str(),fw.c_str(),tag.c_str());
      
  //     const auto _v=
  // 	rawData.find(what);
  //     if(_v==rawData.end())
  // 	{
  // 	  cout<<"List of corr:"<<endl;
  // 	  for(const auto& [key,val] : rawData)
  // 	    cout<<key<<endl;
	  
  // 	  CRASH("Unable to find %s",what.c_str());
  // 	}
      
  //     const auto& v=_v->second;
      
  //     for(size_t iConf=0;iConf<nConfs;iConf++)
  // 	for(size_t t=0;t<T;t++)
  // 	  {
  // 	    union
  // 	    {
  // 	      complex<double> c{};
  // 	      double d[2];
  // 	    };
	    
  // 	    for(size_t ri=0;ri<2;ri++)
  // 	      d[ri]=v[iConf][0][t+T*ri];
	    
  // 	    res[idx({t,iConf})]+=c;
  // 	  }
      
  //     return res;
  //   };
  
  // auto get=
  //   [&idx,
  //    &getRaw,
  //    &nConfs](const auto& rawData,
  // 	      const std::string& tag,
  // 	      const std::string& bw,
  // 	      const std::string& fw,
  // 	      const bool& ri)
  //   {
  //     djvec_t res(T);
      
  //     const auto d=
  // 	getRaw(rawData,tag,bw,fw);
      
  //     for(size_t t=0;t<T;t++)
  // 	{
  // 	  jackknivesFill(nConfs,
  // 			 [&](const size_t& iConf,
  // 			     const size_t& iClust,
  // 			     const double& weight)
  // 			 {
  // 			   const double r=d[idx({t,iConf})].real();
  // 			   const double i=d[idx({t,iConf})].imag();
  // 			   res[t][iClust]+=weight*(ri?i:r);
  // 		       });
  // 	}
  //     res.clusterize(((double)nConfs/njacks));
      
  //     return res;
  //   };
  
  // const auto getMPCAC=
  //   [&rawData,
  //    &get](const std::string& i)
  //   {
  //     const djvec_t AP5P5=get(rawData,"P5P5",i,i,0);
  //     const djvec_t AV0P5=get(rawData,"V0P5",i,i,1);
      
  //     // effective_mass
  //     return (forward_derivative(AV0P5.symmetrized(-1))/AP5P5.symmetrized());
  //   };
  
  // effective_mass(get(rawData,"P5P5","H0","H0",0).symmetrized()).ave_err().write("plots/P5P5_run0.xmg");
  // effective_mass(get(rawData,"P5P5","H1","H1",0).symmetrized()).ave_err().write("plots/P5P5_run1.xmg");
  
  // const djvec_t M0=getMPCAC("H0");
  // M0.ave_err().write("plots/mPCAC_run0.xmg");
  // const djvec_t M1=getMPCAC("H1");
  // M1.ave_err().write("plots/mPCAC_run1.xmg");
  
  // const auto getMPCAC2=
  //   [&rawData,
  //    &get](const std::string& i)
  //   {
  //     const djvec_t AP5P5=get(rawData,"P5P5",i,i,0);
  //     const djvec_t AV0P5=get(rawData,"V0P5",i,i,1);
      
  //     djvec_t out(T);
  //     for(size_t it=0;it<T;it++)
  // 	out[it]=(AV0P5[(it+1)%T]-AV0P5[(it+T-1)%T])/2.0;
  //     // for(size_t it=0;it<T;it++)
  //     // 	out[it]=AV0P5[(it+1)%T]-AV0P5[it];
      
  //     // effective_mass
  //     return out/AP5P5;
  //   };
  
  // getMPCAC2("H0").ave_err().write("plots/mPCAC_Alt.xmg");
  
  
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
