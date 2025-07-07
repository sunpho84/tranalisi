#include "ave_err.hpp"
#include "effective.hpp"
#include "fit.hpp"
#include "functions.hpp"
#include "index.hpp"
#include "jack.hpp"
#include "math.hpp"
#include "meas_vec.hpp"
#include "obs_file.hpp"
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
  djvec_t& A=data.at(#A)[0];		\
  A.symmetrize()
  
  INS(Pi);
  // INS(Pi_dth);
  INS(D);
  INS(D_sme);
  INS(D_sme_sme);
  
#undef INS
  
  // const double th=1e-5;
  // const double kPi=M_PI/L*th;
  // cout<<"k: "<<kPi<<endl;
  
  djack_t ZPiLL,MPi;
  two_pts_fit(ZPiLL,MPi,Pi,T/2,20,30,"plots/Pi.xmg");
  const djack_t ZPiL=sqrt(ZPiLL);
  cout<<"MPi:     "<<smart_print(MPi)<<endl;
  
  const djack_t a=MPi/0.139;
  cout<<"a: "<<a<<endl;
  
  // const djack_t EPiExp=latt_en(MPi,kPi);
  // cout<<"EPi exp: "<<smart_print(EPiExp)<<endl;
  
  // djack_t ZPiLL_dth,EPidTh;
  // two_pts_fit(ZPiLL_dth,EPidTh,Pi_dth,T/2,20,30,"plots/Pi_dth.xmg");
  // cout<<"EPi fit: "<<smart_print(EPidTh)<<endl;
  
  // const djack_t dEPi=EPidTh-MPi;
  // cout<<"Epi-Mpi: "<<smart_print(dEPi)<<", significativity: "<<dEPi.significativity()<<endl;
  
  // (effective_mass(Pi_dth)-effective_mass(Pi)).ave_err().write("plots/pi_th_rat.xmg");
  
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
  // const djack_t Q2=sqr(MD-EPidTh)-sqr(kPi);
  
  cout<<"Q2 max: "<<Q2Max.ave_err()<<endl;
  // cout<<"Q2    : "<<Q2.ave_err()<<endl;
  
  const size_t dT=50;
  
#define INS(A,B,C)							\
  djvec_t A=data.at(#A)[B].subset(C,C+dT)*exp(MPi*C)/ZPiL/ZDS*MD*MPi*4;	\
  for(size_t t=0;t<A.size();t++)					\
    A[t]*=exp(MD*t);							\
  A.ave_err().write("plots/"#A".xmg")
  
  INS(V0_20,0,20);
  INS(V0_30,0,30);
  
  INS(dth_V0_20,1,20);
  INS(V1_20,1,20);
  INS(dth_V1_20,0,20);
  
  INS(dth_V0_30,1,30);
  INS(V1_30,1,30);
  INS(dth_V1_30,0,30);
  
  const djvec_t V0_5_DD=data.at("V0_5_DD")[0];
  const djvec_t V0_10_DD=data.at("V0_10_DD")[0];
  
  const djvec_t V0_20_PP=data.at("V0_20_PP")[0];
  const djvec_t V0_30_PP=data.at("V0_30_PP")[0];
  
  const auto cZ=
    [](const djvec_t& in,
       const djack_t& Z,
       const djack_t& M,
       const size_t& tw)
    {
      djvec_t out(T/2+1-tw);
      for(size_t t=0;t<=T/2-tw;t++)
	out[t]=two_pts_corr_fun(Z,M,T/2.0,t+tw,0)/in[t+tw];
      
      return out;
    };
  
  const djvec_t Z_5_DD=cZ(V0_5_DD,ZDSS,MD,5);
  const djvec_t Z_10_DD=cZ(V0_10_DD,ZDSS,MD,10);
  Z_5_DD.ave_err().write("plots/Z_5_DD.xmg");
  Z_10_DD.ave_err().write("plots/Z_10_DD.xmg");
  
  const djvec_t Z_20_PP=cZ(V0_20_PP,ZPiLL,MPi,20);
  const djvec_t Z_30_PP=cZ(V0_30_PP,ZPiLL,MPi,30);
  Z_20_PP.ave_err().write("plots/Z_20_PP.xmg");
  Z_30_PP.ave_err().write("plots/Z_30_PP.xmg");
  
  const djack_t Z=Z_20_PP[5];
  
  const djvec_t f0_20=V0_20*Z/(MD+MPi);
  const djvec_t f0_30=V0_30*Z/(MD+MPi);
  f0_20.ave_err().write("plots/f0_20.xmg");
  f0_30.ave_err().write("plots/f0_30.xmg");
  
  return 0;
}
