#include <tranalisi.hpp>

using namespace std;

void anSources()
{
  const size_t nSources=2;
  const size_t nConfs=401;
  const size_t T=160;
  const string ens="C80";
  
  // const size_t nSources=8;
  // const size_t nConfs=500;
  // const size_t T=192;
  // const string ens="D96";
  
  index_t id({{"Source",nSources},{"Conf",nConfs},{"T",T}});
  
  vector<double> rawData(id.max());
  
  for(size_t iSource=0;iSource<nSources;iSource++)
    {
      //const string finPath="C80_gen_light_S"+std::to_string(iSource+1)+"_g5.txt";
      const string finPath=ens+"_std_light_S"+std::to_string(iSource+1)+"_g5g4.txt";
      //const string finPath="C80_gen_light_S"+std::to_string(iSource+1)+"_1.txt";
      ifstream fin(finPath);
      
      if(not fin.good())
	CRASH("Error opening %s",finPath.c_str());
      
      for(size_t iConf=0;iConf<nConfs;iConf++)
	{
	  string dum;
	  fin>>dum>>dum;
	  
	  for(size_t t=0;t<T;t++)
	    {
	      double dum;
	      //if(not (fin>>rawData[id({iSource,iConf,t})]>>dum))
	      if(not (fin>>dum>>rawData[id({iSource,iConf,t})]))
		CRASH("Reading conf %zu t %zu",iConf,t);
	    }
	}
    }
  
  vector<double> all(nConfs*T,0.0),unb(nConfs*T,0.0);
  
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      for(size_t t1=0;t1<T;t1++)
	for(size_t dt=0;dt<T;dt++)
	  {
	    const size_t t2=(t1+dt)%T;
	    
	    double s12=0,s1=0,s2=0;
	    for(size_t iSource=0;iSource<nSources;iSource++)
		{
		  const double& x1=rawData[id({iSource,iConf,t1})];
		  const double& x2=rawData[id({iSource,iConf,t2})];
		  s12+=x1*x2;
		  s1+=x1;
		  s2+=x2;
		}
	    
	    all[dt+T*iConf]+=s1*s2;
	    unb[dt+T*iConf]+=s1*s2-s12;
	  }
    }
  all/=-80*80*80*2;
  unb/=-80*80*80*2;
  
  all/=T*sqr(nSources);
  unb/=T*nSources*(nSources-1);
  
  djvec_t a(T,0.0),u(T,0.0);
  
  jackknivesFill(nConfs,
		 [&a,
		  &u,
		  &all,
		  &unb](const size_t& iConf,
			const size_t& iClust,
			const double& w)
		 {
		   for(size_t t=0;t<T;t++)
		     {
		       a[t][iClust]+=all[t+T*iConf]*w;
		       u[t][iClust]+=unb[t+T*iConf]*w;
		     }
		 });
  
  a.clusterize((double)nConfs/njacks).symmetrize();
  u.clusterize((double)nConfs/njacks).symmetrize();
  const djvec_t b=(a-u);
  const djvec_t r=b/u;
  
  a.ave_err().write("plots/all.xmg");
  u.ave_err().write("plots/unb.xmg");
  b.ave_err().write("plots/bias.xmg");
  r.ave_err().write("plots/ratioBias.xmg");
  
  effective_mass(a).ave_err().write("plots/allEff.xmg");
  effective_mass(u).ave_err().write("plots/unbEff.xmg");
  effective_mass(b).ave_err().write("plots/biasEff.xmg");
  
  /////////////////////////////////////////////////////////////////
  
  raw_file_t connF("/home/francesco/QCD/LAVORI/GM3/C.06.80/data/ll_OS_A0A0","r");
  const size_t nC=connF.bin_read<size_t>();
  [[maybe_unused]] const size_t nH=connF.bin_read<size_t>();
  const size_t nT=connF.bin_read<size_t>();

  index_t idC({{"c",nC},{"thP1",nT/2+1}});
  vector<double> rawDataC(idC.max());
  connF.bin_read(rawDataC);
  // if(not connF.feof())
  //   CRASH("feof not reached");
  
  djvec_t conn(T/2+1);
  jackknivesFill(nC,
		 [&conn,
		  &rawDataC,
		  &idC](const size_t& iConf,const size_t& iClust,const double& w)
  {
    for(size_t t=0;t<=T/2;t++)
      conn[t][iClust]+=rawDataC[idC({iConf,t})]*w;
  });
  
  conn.clusterize(double(nC)/njacks);
  conn.ave_err().write("plots/conn.xmg");
  effective_mass(conn).ave_err().write("plots/effConn.xmg");
}

// void anAve()
// {
//   const size_t nConfs=400;
//   const size_t T=160;
  
//   index_t id({{"Conf",nConfs},{"T",T}});
  
//   vector<double> rawData(id.max());
  
//   //const string finPath="C80_std_light_g5g0.txt";
//   const string finPath="C80_gen_light_1.txt";
//   ifstream fin(finPath);
  
//   if(not fin.good())
//     CRASH("Error opening %s",finPath.c_str());
  
//   for(size_t iConf=0;iConf<nConfs;iConf++)
//     {
//       string dum;
//       fin>>dum>>dum;
      
//       for(size_t t=0;t<T;t++)
// 	{
// 	  double dum;
// 	  if(not (fin>>dum>>rawData[id({iConf,t})]))
// 	    CRASH("Reading conf %zu t %zu",iConf,t);
// 	}
//     }
  
//   vector<double> all(nConfs*T,0.0);
  
//   for(size_t iConf=0;iConf<nConfs;iConf++)
//     {
//       for(size_t t1=0;t1<T;t1++)
// 	for(size_t dt=0;dt<T;dt++)
// 	  {
// 	    const size_t t2=(t1+dt)%T;
	    
// 	    const double& x1=rawData[id({iConf,t1})];
// 	    const double& x2=rawData[id({iConf,t2})];
// 	    all[dt+T*iConf]=x1*x2;
// 	  }
//     }
  
//   all/=T;
  
//   djvec_t a(T,0.0);
  
//   jackknivesFill(nConfs,
// 		 [&a,
// 		  &all](const size_t& iConf,
// 			const size_t& iClust,
// 			const double& w)
// 		 {
// 		   for(size_t t=0;t<T;t++)
// 		     {
// 		       a[t][iClust]+=all[t+T*iConf]*w;
// 		     }
// 		 });
  
//   a.clusterize((double)nConfs/njacks);
  
//   a.ave_err().write("plots/allAve.xmg");
  
//   effective_mass(a.symmetrized()).ave_err().write("plots/allAveEff.xmg");
// }

int main()
{
  set_njacks(50);
  
  anSources();
  //anAve();
  
  return 0;
}
