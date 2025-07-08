#include "ave_err.hpp"
#include "effective.hpp"
#include "fit.hpp"
#include "functions.hpp"
#include "index.hpp"
#include "jack.hpp"
#include "math.hpp"
#include "meas_vec.hpp"
#include "obs_file.hpp"
#include "oper.hpp"
#include "raw_file.hpp"
#include <filesystem>
#include <string>
#include <vector>

const size_t L=64;
const size_t T=128;

vector<string> confs;
size_t nConfs;

void setConfs()
{
  for(const auto& entry : filesystem::directory_iterator("data"))
    if(entry.is_directory())
      {
	const string& conf=filesystem::path(entry.path()).filename();
	confs.emplace_back(conf);
      }
  nConfs=confs.size();
  set_njacks(nConfs);
}

int main()
{
  setConfs();
  cout<<"NConfs: "<<nConfs<<endl;
  
  map<string,array<djvec_t,2>> data;
  for(const auto& entry : filesystem::directory_iterator("data/"+confs.front()))
    {
      const string tag=entry.path().filename().string().substr(10);
      vector<vector<double>> rawData;
      for(const filesystem::path conf : confs)
	rawData.push_back(obs_file_t("data"/conf/("mes_contr_"+tag),2,{0,1}).read(T));
      // cout<<tag<<" "<<rawData.back()<<endl;
      
      for(size_t ri=0;ri<2;ri++)
	{
	  djvec_t& d=data[tag][ri];
	  d.resize(T);
	  jackknivesFill(nConfs,[&rawData,
				 &ri,
				 &d](const size_t& iConf,
				     const size_t& iClust,
				     const double& weight)
	  {
	    for(size_t t=0;t<T;t++)
	      d[t][iClust]+=rawData[iConf][t+T*ri]*weight;
	  });
	  d.clusterize((double)nConfs/njacks);
	}
    }
  
  auto get=[&data](const string& name,
		   const bool reim) -> auto&
  {
    auto f=data.find(name);
    if(f==data.end())
      CRASH("unable to find %s",name.c_str());
    
    return f->second[reim];
  };
  
  {
    int i=0;
    for(const auto& [tag,dum] : data)
      {
	if(i++) cout<<",";
	cout<<tag;
      }
    cout<<endl;
  }
  
#define INS(A)				\
  djvec_t& A=get(#A,0);		\
  A.symmetrize()
  
  INS(Pi);
  INS(D);
  INS(D_sme);
  INS(D_sme_sme);
  
#undef INS
  
  djack_t ZPiLL,MPi;
  two_pts_fit(ZPiLL,MPi,Pi,T/2,20,30,"plots/Pi.xmg");
  const djack_t ZPiL=sqrt(ZPiLL);
  cout<<"MPi:     "<<smart_print(MPi)<<endl;
  
  //aggiungere K, Ds
  const djack_t a=MPi/0.139;
  cout<<"a: "<<a<<endl;
  
  djack_t ZDLL,ZDLS,ZDSS,MD;
  two_pts_fit(ZDLL,MD,D,T/2,20,25,"plots/D.xmg");
  two_pts_fit(ZDLS,MD,D_sme,T/2,14,20,"plots/D_sme.xmg");
  two_pts_fit(ZDSS,MD,D_sme_sme,T/2,14,20,"plots/D_sme_sme.xmg");
  const djack_t ZDS=sqrt(ZDSS);
  const djack_t ZDL=sqrt(ZDLL);
  cout<<"ZDL: "<<ZDL.ave_err()<<endl;
  cout<<"ZDS: "<<ZDS.ave_err()<<endl;
  cout<<"ZDLS: "<<ZDLS.ave_err()<<" ZDL*ZDS: "<<(ZDL*ZDS).ave_err()<<endl;
  cout<<"MD: "<<MD.ave_err()<<endl;
  
  const djack_t Q2Max=sqr(MD-MPi);
  
  cout<<"Q2 max: "<<Q2Max.ave_err()<<endl;
  
  const size_t dT=50;
  
#define INS(A,B,C)							\
  djvec_t A=get(#A,B).subset(C,C+dT)*exp(MPi*C)/ZPiL/ZDS*MD*MPi*4;	\
  for(size_t t=0;t<A.size();t++)					\
    A[t]*=exp(MD*t);							\
  A.ave_err().write("plots/"#A".xmg")
  
  INS(V0_30,0,30);
  INS(V0_40,0,40);
  
  INS(S0_30,0,30);
  INS(S0_40,0,40);
  
  // const djvec_t V0_5_DD=get("V0_5_DD",0);
  // const djvec_t V0_10_DD=get("V0_10_DD",0);
  
  // const djvec_t V0_30_PP=get("V0_30_PP",0);
  // const djvec_t V0_40_PP=get("V0_40_PP",0);
  
  // const auto cZ=
  //   [](const djvec_t& in,
  //      const djack_t& Z2,
  //      const djack_t& M,
  //      const size_t& tw)
  //   {
  //     djvec_t out(T/2+1-tw);
  //     for(size_t t=0;t<=T/2-tw;t++)
  // 	out[t]=Z2/(2*M)*exp(-M*(t+tw))/in[t+tw];
      
  //     return out;
  //   };
  
  // const djvec_t Z_5_DD=cZ(V0_5_DD,ZDSS,MD,5);
  // const djvec_t Z_10_DD=cZ(V0_10_DD,ZDSS,MD,10);
  // Z_5_DD.ave_err().write("plots/Z_5_DD.xmg");
  // Z_10_DD.ave_err().write("plots/Z_10_DD.xmg");
  
  // const djvec_t Z_30_PP=cZ(V0_30_PP,ZPiLL,MPi,30);
  // const djvec_t Z_40_PP=cZ(V0_40_PP,ZPiLL,MPi,40);
  // Z_30_PP.ave_err().write("plots/Z_30_PP.xmg");
  // Z_40_PP.ave_err().write("plots/Z_40_PP.xmg");
  
  // const djack_t Z=Z_30_PP[5];
  // cout<<"Z: "<<smart_print(Z)<<endl;
  
  const djvec_t f0_30=S0_30*(0.231567-0.00072)/(sqr(MD)-sqr(MPi));
  const djvec_t f0_40=S0_40*(0.231567-0.00072)/(sqr(MD)-sqr(MPi));
  f0_30.ave_err().write("plots/f0_30.xmg");
  f0_40.ave_err().write("plots/f0_40.xmg");
  
  const djvec_t Z_30=S0_30/V0_30*(0.231567-0.00072)/(MD-MPi);
  const djvec_t Z_40=S0_40/V0_40*(0.231567-0.00072)/(MD-MPi);
  Z_30.ave_err().write("plots/Z_30.xmg");
  Z_40.ave_err().write("plots/Z_40.xmg");
  
  INS(dth_V1_30,0,30);
  INS(dth_V1_40,0,40);
  
  const djvec_t fp_30_Pi_D=(V0_30+(MD-MPi)*dth_V1_30)*Z_30/(2*MD);
  const djvec_t fp_40_Pi_D=(V0_40+(MD-MPi)*dth_V1_40)*Z_40/(2*MD);
  fp_30_Pi_D.ave_err().write("plots/fp_30_Pi_D.xmg");
  fp_40_Pi_D.ave_err().write("plots/fp_40_Pi_D.xmg");
  
  INS(V1_dth_30,0,30);
  INS(V1_dth_40,0,40);
  
  const djvec_t fp_30_D_Pi=(V0_30+(MPi-MD)*V1_dth_30)*Z_30/(2*MPi);
  const djvec_t fp_40_D_Pi=(V0_40+(MPi-MD)*V1_dth_40)*Z_40/(2*MPi);
  fp_30_D_Pi.ave_err().write("plots/fp_30_D_Pi.xmg");
  fp_40_D_Pi.ave_err().write("plots/fp_40_D_Pi.xmg");
  
  return 0;
}
