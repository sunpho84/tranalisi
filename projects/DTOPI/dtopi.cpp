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
#include "random.hpp"
#include "raw_file.hpp"
#include <filesystem>
#include <string>
#include <vector>

const size_t L=64;
const size_t T=128;

const int nls=1;
const size_t nc=1;
const double ml=0.00066690;
const double ms=0.018267;


vector<string> confs;
size_t nConfs;

void setConfs()
{
  for(const auto& entry : filesystem::directory_iterator("data"))
    if(entry.is_directory())
      //if(filesystem::exists(entry.path()/"mes_contr_3pUnsme"))
      {
	const string& conf=filesystem::path(entry.path()).filename();
	confs.emplace_back(conf);
      }
  nConfs=confs.size();
  set_njacks(70);
}

int main()
{
  setConfs();
  cout<<"NConfs: "<<nConfs<<endl;
  
  map<string,array<djvec_t,2>> data;
  for(const char* entry : {"2p","3p"// ,"3pUnsme"
    })
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
  
  const map<string,ave_err_t> liter=
    {{"fnalpD0ToPi",{3.119,0.057}},
     {"fnal0D0ToPi",{1.2783,0.0061}},
     {"etmpD0ToPi",{2.130,0.096}},
     {"etm0D0ToPi",{1.134,0.049}},
     {"fnalpD0ToK",{1.451,0.017}},
     {"fnal0D0ToK",{1.0240,0.0021}},
     {"etmpD0ToK",{1.336,0.054}},
     {"etm0D0ToK",{0.979,0.019}},
     {"fnalpD0sToK",{1.576,0.013}},
     {"fnal0D0sToK",{0.9843,0.0018}}};
  
  auto get=[&data](const string& name,
		   const bool reim) -> auto&
  {
    auto f=data.find(name);
    if(f==data.end())
      CRASH("unable to find %s",name.c_str());
    
    return f->second[reim];
  };
  
  djack_t ZPLList[2],MPList[2];
  for(size_t ils=0;ils<nls;ils++)
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
  if(not file_exists("input"))
    CRASH("input file not present");
  input_file_t input("input");
  const size_t tw=input.read<size_t>("tw");
  const vector<size_t> twList={tw};
  const size_t ntw=twList.size();
  const size_t nsme=1;
  
  const index_t iDs({{"ls",nls},{"c",nc},{"sme",nsme}});
  
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
  for(size_t ils=0;ils<nls;ils++)
    for(size_t ic=0;ic<nc;ic++)
      {
	const string ls=((ils==0)?"l":"s");
	
	djack_t& M=MD(ils,ic,0);
	const djvec_t D=get("S_m"+ls+"__S_mc"+to_string(ic),0).symmetrized();
	const string Dname="D"+to_string(ic)+((ils==1)?"s":"");
	
	djack_t ZDLL;
	two_pts_fit(ZDLL,M,D,T/2,27,30,"plots/"+Dname+".xmg");
	
	for(size_t is=0;is<nsme;is++)
	  {
	    const djvec_t D_sme=get("S_m"+ls+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	    const djvec_t D_sme_sme=get("S_m"+ls+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	    
	    djack_t ZDLS,ZDSS;
	    two_pts_fit(ZDLS,M,D_sme,T/2,22,28,"plots/"+Dname+"_sme"+to_string(is)+".xmg");
	    two_pts_fit(ZDSS,M,D_sme_sme,T/2,22,28,"plots/"+Dname+"_sme"+to_string(is)+"_sme"+to_string(is)+".xmg");
	    
	    const djack_t& l=ZDS(ils,ic,is)=sqrt(ZDSS);
	    const djack_t& s=ZDL(ils,ic,is)=sqrt(ZDLL);
	    
	    cout<<"M"<<Dname<<": "<<(M/a).ave_err()<<endl;
	    cout<<"ZDL: "<<l.ave_err()<<endl;
	    cout<<"ZDS: "<<s.ave_err()<<endl;
	    cout<<"ZDLS: "<<ZDLS.ave_err()<<" ZDL*ZDS: "<<(l*s).ave_err()<<endl;
	    
	    // const size_t t0=5;
	    // auto [eig,recastEigvec,origEigvec]=
	    //   gevp({D,D_sme,D_sme,D_sme_sme},t0);
	    
	    // effective_mass(eig[0]).ave_err().write("plots/"+Dname+"_eig1.xmg");
	    // effective_mass(eig[1]).ave_err().write("plots/"+Dname+"_eig2.xmg");
	  }
	
      }
  
  if(ntw>1 or nsme>1)
    CRASH("Too many tw or sme");

  std::vector<pair<size_t,size_t>> lsCombo{{0,0}};
  if(nls>1)
    for(const auto& [l,s] : {std::make_pair(0,1),{1,0}})
      lsCombo.emplace_back(l,s);
  
  for(const auto& ilsSS : lsCombo)
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
      
      const size_t dT=60;
      
      const size_t tmin3=27,tmax3=35;
      djvec_t f0(nc),fp(nc),ft(nc);
      for(size_t ic=0;ic<nc;ic++)
	{
	  const string Dname="D"+to_string(ic)+((ilsSpe==1)?"s":"");
	  const double mc=mcList[ic];
	  
	  djack_t& M=MD(ilsSpe,ic,0);
	  const djack_t Q2Max=sqr(M-MP)/sqr(a);
	  cout<<Dname<<" -> "<<Pname<<" Q2 max: "<<Q2Max.ave_err()<<endl;
	  
	  for(size_t is=0;is<nsme;is++)
	    {
	      const djvec_t D_sme=get("S_m"+lsSpe+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	      const djvec_t D_sme_sme=get("S_m"+lsSpe+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	      
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
		  const djvec_t cS0=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  
		  const string suff="_"+Dname+"_sme"+to_string(is)+"_tw"+to_string(tw)+"to"+Pname+".xmg";
		  
		  // const djack_t cc=ZDS(ilsSpe,ic,is)/ZDL(ilsSpe,ic,is);
		  // const djvec_t cV0Unsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_V0_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  // const djvec_t cdthVUnsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_dth_V1_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  // const djvec_t cdthTUnsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_dth_T1_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  // const djvec_t cS0Unsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  
		  const djvec_t cf0=cS0*(mc-mseq)/(sqr(M)-sqr(MP));
		  f0[ic]=constant_fit(cf0,tmin3,tmax3);
		  forward_derivative(-log(cf0)).ave_err().write("plots/effMassF0"+suff);
		  
		  const djvec_t cZ=cS0/cV0*(mc-mseq)/(M-MP)*0+0.7;
		  cZ.ave_err().write("plots/Z"+suff);
		  
		  const djvec_t cfp=(cV0+(M-MP)*cdthV)*cZ/(2*M);
		  fp[ic]=constant_fit(cfp,tmin3,tmax3,"plots/fp"+suff);
		  // const djvec_t cfpUnsme=(cV0Unsme+(M-MP)*cdthVUnsme)*cZ/(2*M);
		  
		  // grace_file_t fpCompa("plots/fpCompa"+suff);
		  // fpCompa.write_vec_ave_err(cfp.ave_err());
		  // fpCompa.set_legend("Smeared");
		  // fpCompa.set_all_colors(grace::BLUE);
		  // fpCompa.write_vec_ave_err(cfpUnsme.ave_err());
		  // fpCompa.set_legend("Local");
		  
		  jack_fit_t fitter;
		  djvec_t pars(3);
		  fitter.add_fit_par(pars[0],"C",4,0.1);
		  fitter.add_fit_par(pars[1],"dM",MPList[0].ave(),0.1);
		  fitter.add_fit_par(pars[2],"dZ",-2.6,0.1);
		  
		  const size_t t2s=15;
		  for(size_t it=t2s;it<=60;it++)
		    fitter.add_point(cfp[it],[t=(double)it](const vector<double>& pars,
							    const size_t iel)
		    {
		      return pars[0]+pars[2]*exp(-pars[1]*t);
		    });
		  fitter.fit();
		  for(const auto& p : pars)
		    cout<<""<<p.ave_err()<<endl;
		  
		  cout<<"dM: "<<(pars[1]/a).ave_err()<<endl;
		  
		  grace_file_t plot2sFitFp("plots/2sfp"+suff);
		  plot2sFitFp.write_vec_ave_err(cfp.ave_err());
		  for(const auto& [b,e] : {std::make_pair(0,t2s),{t2s,60}})
		    plot2sFitFp.write_polygon([&pars](const double& t)
		    {
		      return pars[0]+pars[2]*exp(-pars[1]*t);
		    },b,e);
		  
		  grace_file_t plot2sFit2subFp("plots/2sfp2sub"+suff);
		  djvec_t i=cfp;
		  for(size_t t=0;t<cfp.size();t++)
		    i[t]-=pars[2]*exp(-pars[1]*t);
		  plot2sFit2subFp.write_vec_ave_err(i.ave_err());
		  plot2sFit2subFp.write_polygon([&pars](const double& t)
		    {
		      return pars[0];
		    },0,60);
		  
		  cout<<" t dom: "<<(-log(abs(pars[0]/pars[2]/100))/pars[1]).ave_err()<<endl;
		  
		  const djvec_t cfpfr0=(cV0+(M-MP)*cdthV)*cZ/(2*M)/cf0;
		  constant_fit(cfpfr0,tmin3,tmax3,"plots/fpfr0"+suff);
		  
		  djvec_t cft;
		  
		  // const djvec_t cdthT=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_dth_T1_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  // cft=cdthT*(sqr(M)+sqr(MP))/(2*M);
		  // ft[ic]=constant_fit(cft,tmin3,tmax3,"plots/ft"+suff);
		  
		  for(const auto&[tag,c,f] :
			{make_tuple("p",cfp,fp[ic]),
			 /*      */{"0",cf0,f0[ic]},
			 /*      */{"t",cft,ft[ic]}
			})
		    {
		      grace_file_t plotf("plots/f"+(tag+suff));
		      plotf.set_title((string)tag+"("+Dname+"->"+Pname+")");
		      plotf.set_xaxis_label("t/a");
		      plotf.write_vec_ave_err(c.ave_err());
		      plotf.write_constant_band(tmin3,tmax3,f);
		      for(const string coll : {"etm","fnal"})
			if(const auto i=liter.find(coll+tag+Dname+"To"+Pname);i!=liter.end())
			  {
			    plotf.set_legend(coll);
			    plotf.write_constant_band(0,dT,djack_t(gauss_filler_t(i->second,12432)));
			  }
			// else
			//   cout<<coll+tag+Dname+"To"+Pname<<" not found in the literature"<<endl;
		    }
		  
		  cfp.bin_write("cfp.dat");
		  cf0.bin_write("cf0.dat");
		}
	    }
	}
      
      const string Dname=
	(string)"D"+((ilsSpe==1)?"s":"");
      auto plot=
	[&](const string& tag,
	      const djvec_t& f)
	{
	  const string decay=tag+Dname+"to"+Pname;
	  grace_file_t plotf("plots/f"+decay+".xmg");
	  for(size_t ic=0;ic<nc;ic++)
	    plotf.write_ave_err(MD(0,ic,0).ave(),f[ic].ave_err());
	};
      
      for(const auto&[tag,f] : {make_pair("p",fp),{"0",f0},{"t",ft}})
	plot(tag,f);
    }
  
  return 0;
}
