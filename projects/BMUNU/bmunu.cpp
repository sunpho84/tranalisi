#include <tranalisi.hpp>

int main()
{
  set_njacks(20);
  
  const size_t T=128;
  const size_t nMes2pt=3;
  const size_t nSme=2;
  const size_t nConfs=200;
  const size_t nMb=5;
  
  const size_t nTh=5;
  const size_t nQIns=2;
  const size_t nVmu=2;
  const size_t nFB=2;
  
  const index_t id2({{"Mb",nMb},{"Conf",nConfs},{"Sme",nSme},{"Mes",nMes2pt},{"t",T}});
  cout<<id2.max()*sizeof(double)<<endl;
  
  vector<double> mes2ptData(id2.max());
  raw_file_t mes2ptFile("data/mes_2pt.dat","r");
  mes2ptFile.bin_read(mes2ptData);
  
  auto load2Pts=
    [&](const size_t& iMb,
	const size_t& iSme,
	const size_t& iMes2pt)
    {
      djvec_t corr(T);
      
      jackknivesFill(nConfs,[&](const size_t& iConf,
				const size_t& iClust,
				const double& w)
      {
	for(size_t t=0;t<T;t++)
	  corr[t][iClust]+=w*mes2ptData[id2({iMb,iConf,iSme,iMes2pt,t})];
      });
      
      corr.clusterize((double)nConfs/njacks);
      
      return corr;
    };
  
  const djvec_t pionCorr=
    load2Pts(0,0,0).symmetrized();
  
  effective_mass(pionCorr).ave_err().write("plots/pion.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  const index_t id3({{"Mb",nMb},{"Conf",nConfs},{"Th",nTh},{"QIns",nQIns},{"FB",nFB},{"Vmu",nVmu},{"t",T},{"RI",2}});
  cout<<id3.max()*sizeof(double)<<endl;
  
  vector<double> mes3ptData(id3.max());
  raw_file_t mes3ptFile("data/mes_3pt_A1P5.dat","r");
  mes3ptFile.bin_read(mes3ptData);
  
  auto load3Pts=
    [&](const size_t& iMb,
	const size_t& iTh,
	const size_t& iQIns,
	const size_t& fb,
	const size_t& mu,
	const size_t& ri)
    {
      const double charge=-1.0/3;
      djvec_t corr(T);
      
      jackknivesFill(nConfs,[&](const size_t& iConf,
				const size_t& iClust,
				const double& w)
      {
	for(size_t t=0;t<T;t++)
	  corr[t][iClust]+=w*mes3ptData[id3({iMb,iConf,iTh,iQIns,fb,mu,t,ri,})];
      });
      
      corr.clusterize((double)nConfs/njacks);
      
      return corr*charge;
    };
  
  enum{LIGHT_INS,HEAVY_INS};
  enum{FORW,BACK};
  enum{V1,V2};
  
  const size_t iMb=1,iTh=1,iQins=LIGHT_INS,fb=FORW,mu=V1,ri=RE;
  const djvec_t A1P5=
    load3Pts(iMb,iTh,iQins,fb,mu,ri);
  
  A1P5.ave_err().write("plots/A1P5.xmg");
  
  return 0;
}
