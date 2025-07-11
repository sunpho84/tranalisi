#include "ave_err.hpp"
#include "effective.hpp"
#include "fit.hpp"
#include "functions.hpp"
#include "grace.hpp"
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

const double ml=0.00066690;
const double ms=0.018267;

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
  for(const char* entry : {"2p","3p"})
    {
      map<string,vector<vector<double>>> rawData;
      
      for(const filesystem::path conf : confs)
	{
	  raw_file_t file("data"/conf/("mes_contr_"+(string)entry),"r");
	  char line[1024];
	  auto readLine=[&file,
			 &line]()
	  {
	    bool r=false;
	    while((not r) and (not file.feof()))
	      {
		file.get_line(line);
		for(char* c=line;*c!='\0' and not r;c++)
		  r|=(*c!=' ' and *c!='\0');
	      };
	    
	    return r;
	  };
	  
	  char a[100],b[100];
	  while(readLine())
	    {
	      if(sscanf(line," # Contraction of %s ^ \\dag and %s",a,b)!=2)
		CRASH("Unable to parse %s",line);
	      
	      const string tag=(string)a+"__"+b;
	      
	      if(not readLine())
		CRASH("Unable to get the corr name");
	      
	      vector<double>& data=rawData[tag].emplace_back(2*T);
	      for(size_t t=0;t<T;t++)
		{
		  if(not readLine())
		    CRASH("Unable to read time %zu for contr %s %s",t,a,b);
		  
		  double r,i;
		  if(sscanf(line,"%lg %lg",&r,&i)!=2)
		    CRASH("Unable to convert %s to two doubles",line);
		  
		  data[t]=r;
		  data[t+T]=i;
		}
	      // cout<<tag<<" "<<rawData.back()<<endl;
	    }
	}
      
      for(const auto& it : rawData)
	for(size_t ri=0;ri<2;ri++)
	  {
	    const string& tag=it.first;
	    const vector<vector<double>>& v=it.second;
	    djvec_t& d=data[tag][ri];
	    d.resize(T);
	    jackknivesFill(nConfs,[&v,
				   &ri,
				   &d](const size_t& iConf,
				       const size_t& iClust,
				       const double& weight)
	    {
	      for(size_t t=0;t<T;t++)
		d[t][iClust]+=v[iConf][t+T*ri]*weight;
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
  
  djack_t ZPLList[2],MPList[2];
  for(size_t ils=0;ils<2;ils++)
    {
      const string ls=(ils==0)?"l":"s";
      const djvec_t P=get("S_m"+ls+"__S_ml",0).symmetrized();
      djack_t ZPLL;
      const string Pname=ils?"K":"Pi";
      two_pts_fit(ZPLL,MPList[ils],P,T/2,20,30,"plots/"+Pname+".xmg");
      ZPLList[ils]=sqrt(ZPLL);
      cout<<"MP"+ls+":     "<<smart_print(MPList[ils])<<endl;
    }
  
  const djack_t fP=ZPLList[0]*2*ml/sqr(MPList[0]);
  
  const djack_t a=MPList[0]/0.135;
  cout<<"a: "<<a<<endl;
  cout<<"fpi: "<<(fP/a).ave_err()<<" GeV"<<endl;
  
  const djack_t xi=MPList[0]/fP;
  cout<<"xi: "<<xi.ave_err()<<endl;
  
  const vector<double> mcList={0.231567,0.287798,0.357684,0.444541,0.5524889};
  const vector<size_t> twList={30};//{30,40};
  const size_t ntw=twList.size();
  const size_t nsme=1;
  const size_t nc=mcList.size();
  
  const index_t iDs({{"ls",2},{"c",nc},{"sme",nsme}});
  
  auto getObs=
    [&iDs]()
    {
      return
	[iDs,
	 data=djvec_t(iDs.max())](const size_t& ils,
				  const size_t& ic,
				  const size_t& isme) mutable -> djack_t&
	{
	  return data[iDs({ils,ic,isme})];
	};
    };
  auto MD=getObs();
  auto ZDL=getObs();
  auto ZDS=getObs();
  for(size_t ils=0;ils<2;ils++)
    for(size_t ic=0;ic<nc;ic++)
      {
	const string ls=((ils==0)?"l":"s");
	
	djack_t& M=MD(ils,ic,0);
	const djvec_t D=get("S_m"+ls+"__S_mc"+to_string(ic),0).symmetrized();
	const string Dname="D"+to_string(ic)+((ils==1)?"s":"");
	
	djack_t ZDLL;
	two_pts_fit(ZDLL,M,D,T/2,20,25,"plots/"+Dname+".xmg");
	cout<<"M"<<Dname<<": "<<(M/a).ave_err()<<endl;
	
	for(size_t is=0;is<nsme;is++)
	  {
	    const djvec_t D_sme=get("S_m"+ls+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	    const djvec_t D_sme_sme=get("S_m"+ls+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	    
	    djack_t ZDLS,ZDSS;
	    two_pts_fit(ZDLS,M,D_sme,T/2,14,20,"plots/"+Dname+"_sme"+to_string(is)+".xmg");
	    two_pts_fit(ZDSS,M,D_sme_sme,T/2,14,20,"plots/"+Dname+"_sme"+to_string(is)+"_sme"+to_string(is)+".xmg");
	    
	    djack_t& l=ZDS(ils,ic,is)=sqrt(ZDSS);
	    djack_t s=ZDL(ils,ic,is)=sqrt(ZDLL);
	    
	    cout<<"ZDL: "<<l.ave_err()<<endl;
	    cout<<"ZDS: "<<s.ave_err()<<endl;
	    cout<<"ZDLS: "<<ZDLS.ave_err()<<" ZDL*ZDS: "<<(l*s).ave_err()<<endl;
	    
	  }
      }
  
  if(ntw>1 or nsme>1)
    CRASH("Too many tw or sme");
  
  for(const auto& ilsSS : {pair<size_t,size_t>{0,0},{0,1},{1,0}})
    {
      const size_t& ilsSeq=ilsSS.first;
      const size_t& ilsSpe=ilsSS.second;
      
      const string lsSeq=((ilsSeq==0)?"l":"s");
      const string lsSpe=((ilsSpe==0)?"l":"s");
      
      const size_t i=ilsSeq+ilsSpe;
      const djack_t MP=MPList[i];
      const djack_t ZP=ZPLList[ilsSeq+ilsSpe];
      const string Pname=(ilsSpe==1 or ilsSeq==1)?"K":"Pi";
      
      const double mseq=(ilsSeq==0)?ml:ms;
      
      djvec_t f0(nc),fp(nc),ft(nc);
      for(size_t ic=0;ic<nc;ic++)
	{
	  const string Dname="D"+to_string(ic)+((ilsSpe==1)?"s":"");
	  const double mc=mcList[ic];
	  
	  djack_t& M=MD(ilsSpe,ic,0);
	  const djack_t Q2Max=sqr(M-MP);
	  cout<<Dname<<" -> "<<Pname<<" Q2 max: "<<Q2Max.ave_err()<<endl;
	  
	  for(size_t is=0;is<nsme;is++)
	    {
	      const djvec_t D_sme=get("S_m"+lsSpe+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	      const djvec_t D_sme_sme=get("S_m"+lsSpe+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	      
	      const size_t dT=50;
	      
	      for(size_t itw=0;itw<ntw;itw++)
		{
		  const size_t tw=twList[itw];
		  auto g=
		    [&](const string& name,
			const bool ri)
		    {
		      djvec_t A=get(name,ri).subset(tw,tw+dT)*exp(MP*tw)/ZP/ZDS(ilsSpe,ic,is)*M*MP*4;
		      for(size_t t=0;t<A.size();t++)
			A[t]*=exp(M*t);
		      A.ave_err().write("plots/"+name+".xmg");
		      
		      return A;
		    };
		  
		  const djvec_t cV0=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_V0_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  const djvec_t cdthV=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_dth_V1_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  const djvec_t cdthT=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_dth_T1_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  const djvec_t cS0=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  
		  const string suff="_"+Dname+"_sme"+to_string(is)+"_tw"+to_string(tw)+"to"+Pname+".xmg";
		  
		  const djvec_t cf0=cS0*(mc-mseq)/(sqr(M)-sqr(MP));
		  f0[ic]=constant_fit(cf0,20,25,"plots/f0"+suff);
		  
		  const djvec_t cZ=cS0/cV0*(mc-mseq)/(M-MP);
		  cZ.ave_err().write("plots/Z"+suff);
		  
		  const djvec_t cfp=(cV0+(M-MP)*cdthV)*cZ/(2*M);
		  fp[ic]=constant_fit(cfp,20,25,"plots/fp"+suff);
		  
		  const djvec_t cft=cdthT*(sqr(M)+sqr(MP))/(2*M);
		  ft[ic]=constant_fit(cft,20,25,"plots/ft"+suff);
		}
	    }
	}
      
      const string Dname=(string)"D"+((ilsSpe==1)?"s":"");
      auto plot=
	[&Dname,
	 &Pname,
	 &nc,
	 &MD](const string& tag,
	      const djvec_t& f)
	{
	  grace_file_t plotf("plots/f"+tag+Dname+"to"+Pname+".xmg");
	  for(size_t ic=0;ic<nc;ic++)
	    plotf.write_ave_err(MD(1,ic,0).ave(),f[ic].ave_err());
	};
      
      for(const auto&[tag,f] : {make_pair("p",fp),{"0",f0},{"t",ft}})
	plot(tag,f);
    }
  
  return 0;
}
