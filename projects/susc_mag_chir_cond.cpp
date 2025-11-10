#include <tranalisi.hpp>


int main()
{
  const size_t T=128;
  
  ifstream ifiDisc("/home/francesco/QCD/LAVORI/SUSC_MAG_ORSAY/new/std_g5.txt");
  
  const size_t nConfsDisco=749;
  map<string,double> disc;
  for(size_t iconf=0;iconf<nConfsDisco;iconf++)
    {
      string dum,conf;
      ifiDisc>>dum>>conf;
      
      double s=0;
      for(size_t t=0;t<T;t++)
	{
	  double r,i;
	  if(not (ifiDisc>>r>>i))
	    CRASH("aaah");
	  s+=i/2;
	}
      disc[conf]=s;
    }
  
  /////////////////////////////////////////////////////////////////
  
  ifstream ifiConnConfsList("/home/francesco/QCD/LAVORI/SUSC_MAG_ORSAY/new/confsListLL");
  string a;
  vector<string> connConfsList;
  while(ifiConnConfsList>>a)
    connConfsList.emplace_back(a);
  cout<<"ConfsListsize: "<<connConfsList.size()<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  raw_file_t ifiConn("/home/francesco/QCD/LAVORI/GM3/B.72.64/data/ll_TM_VKTK","r");
  const size_t nConnConfs=ifiConn.bin_read<size_t>();
  
  ifiConn.bin_read<size_t>();
  ifiConn.bin_read<size_t>();
  
  const index_t iRaw({{"nConfs",nConnConfs},{"T/2+1",T/2+1}});
  vector<double> conn(iRaw.max());
  ifiConn.bin_read(conn);
  
  vector<pair<size_t,string>> refinedList;
  for(size_t iConf=0;iConf<nConnConfs;iConf++)
    if(auto conf=connConfsList[iConf];disc.find(conf)!=disc.end())
      refinedList.emplace_back(iConf,conf);
  
  // cout<<"Refined list: "<<refinedList.size()<<endl;
  // for(const auto& c : refinedList)
  //   cout<<c.second<<endl;
  
  /////////////////////////////////////////////////////////////////
  
  const size_t nConfs=refinedList.size();
  njacks=50;
  
  djvec_t corr(T/2+1);
  djvec_t corrS(T/2+1);
  djack_t S;
  
  jackknivesFill(nConfs,[&](const size_t& jConf,
			    const size_t& iClust,
			    const double& w)
  {
    const auto& [iConf,name]=refinedList[jConf];
    const double& d=disc[name];
    
    S[iClust]+=w*d;
    for(size_t t=0;t<=T/2;t++)
      {
	corr[t][iClust]+=w*conn[iRaw({iConf,t})];
	corrS[t][iClust]+=w*conn[iRaw({iConf,t})]*d;
      }
  });
  
  const double clustSize=nConfs/(double)njacks;
  S.clusterize(clustSize);
  corr.clusterize(clustSize);
  corrS.clusterize(clustSize);
  
  cout<<S.ave_err()<<endl;
  corrS.ave_err().write("corrSuns.xmg");
  corrS-=corr*S;
  corr.ave_err().write("corr.xmg");
  corrS.ave_err().write("corrS.xmg");
  
  // {
  //   djvec_t res(T/2+1);
  //   jackknivesFill(nConfs,
  // 		   [&](const size_t& iConf,
  // 		       const size_t& iJack,
  // 		       const double& w)
  // 		   {
  // 		     for(size_t t=0;t<=T/2;t++)
  // 		       res[t][iJack]+=w*conn[iRaw({iConf,t})];
  // 		   });
  //   res.clusterize((double)nConfs/njacks);
  // }

  // djvec_load
  
  return 0;
}
