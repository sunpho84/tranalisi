#include <tranalisi.hpp>

size_t T,L;

djvec_t load(const string name,const size_t id)
{
  return read_djvec("jacks/"+name,T,id);
}

size_t tMinPi,tMaxPi;
size_t tMinK,tMaxK;

int main(int narg,char **arg)
{
  raw_file_t input("input","r");
  T=input.read<size_t>("T");
  L=T/2;
  set_njacks(input.read<size_t>("NJacks"));
  
  //size_t tMinPi=8,tMaxPi=17;
  //size_t tMinK=9,tMaxK=17;
  
  const size_t tMinPi=input.read<size_t>("TIntPi");
  const size_t tMaxPi=input.read<size_t>();
  const size_t tMinK=input.read<size_t>("TIntK");
  const size_t tMaxK=input.read<size_t>();
  
  [[ maybe_unused ]]
  const double a=0.4568;
  
  const djvec_t PiLoc=load("S0_M0_R0_0_S0_M0_R0_0_P5P5",0).symmetrized();
  const djvec_t PiSme=load("S0_M0_R0_0_SM_S0_M0_R0_0_SM_P5P5",0).symmetrized();
  djack_t ZPL,ZPS,MP;
  two_pts_SL_fit(ZPS,ZPL,MP,PiLoc,PiSme,T/2,tMinPi,tMaxPi,"plots/MP.xmg");
  cout<<"aMP: "<<smart_print(MP)<<endl;
  // cout<<"ZPL: "<<smart_print(ZPL)<<endl;
  // cout<<"ZPS: "<<smart_print(ZPS)<<endl;
  const djack_t aFp=2*0.0040*ZPL/MP/MP;
  const djack_t fP=aFp/a;
  cout<<"fP: "<<smart_print(fP)<<endl;
  // cout<<"mP: "<<smart_print((djack_t)(MP/a))<<endl;
  
  const double piMin=2*M_PI/L;
  cout<<"Pi min: "<<piMin<<endl;
  const djack_t EpiMin=sqrt(MP*MP+piMin*piMin);
  cout<<"E pi in motion: "<<smart_print(EpiMin)<<endl;
  
  const djvec_t KLoc=load("S0_M1_R0_0_S0_M0_R0_0_P5P5",0).symmetrized();
  const djvec_t KSme=load("S0_M1_R0_0_SM_S0_M0_R0_0_SM_P5P5",0).symmetrized();
  djack_t ZKL,ZKS,MK;
  two_pts_SL_fit(ZKS,ZKL,MK,KLoc,KSme,T/2,tMinK,tMaxK,"plots/MK.xmg");
  cout<<"aMK: "<<smart_print(MK)<<endl;
  
  const size_t tSepList[3]={20,24,28};
  vector<djvec_t> KVP_corr(3,djvec_t(T));
  vector<djvec_t> dt(3,djvec_t(T));
  djvec_t Mel(3);
  
  for(size_t iTsep=0;iTsep<3;iTsep++)
    {
      const size_t tSep=tSepList[iTsep];
      
      for(size_t it=0;it<=tSep;it++)
	{
	  const double t=it;
	  dt[iTsep][it]=-ZPS*ZKS*exp(-MK*t)*exp(-MP*(tSep-t))/(4*MK*MP);
	  //dt[it]-=Z2P*Z2K*exp(-MP*(T-tSep))*exp(-MK*t);
	}
      for(size_t it=tSep;it<T;it++)
	{
	  const double t=it;
	  dt[iTsep][it]=ZPS*ZKS*exp(-MK*(T-t))*exp(-MP*(t-tSep))/(4*MK*MP);
	}
      
      // const djvec_t PVP=load(combine("S0_M0_R0_0_S1_M0_R0_%zu_V0P5",tSep),0)/dt;
      // PVP.ave_err().write(combine("plots/PVP_%zu.xmg",tSep));
      
      KVP_corr[iTsep]=load(combine("S0_M1_R0_0_S1_M0_R0_%zu_V0P5",tSep),0);
      const djvec_t KVP=KVP_corr[iTsep]/dt[iTsep];
      Mel[iTsep]=constant_fit(KVP,tSep/2-2,tSep/2+2,combine("plots/KVP_%zu.xmg",tSep));
      
      cout<<"Mel: "<<Mel[iTsep]<<endl;
      
      const djvec_t KVP_exc=KVP_corr[iTsep]-Mel[iTsep]*dt[iTsep];
      const djvec_t KVP_exc_b=log(KVP_exc/KVP_exc.shift(1));
      
      jack_fit_t fit;
      
      djvec_t pFit(3);
      fit.add_fit_par(pFit[0],"A0",-1e-6,0.1e-6);
      fit.add_fit_par(pFit[1],"A1",+1e-6,0.1e-6);
      fit.add_fit_par(pFit[2],"MK",+0.3,0.1);

      const size_t dt=1;
      for(size_t it=tSep+dt;it<T-dt;it++)
	{
	  fit.add_point(KVP_exc[it],
			[=]
			(const vector<double>& p,int ijack)
			{
			  return p[0]*exp(-p[2]*(it-tSep))+p[1]*exp(-p[2]*(T-it));
			});
	}
      
      auto status=fit.fit();
      
      cout<<pFit[2].ave_err()<<endl;
      
      grace_file_t KVP_exc_plot(combine("plots/KVP_exc_%zu.xmg",tSep));
      KVP_exc_plot.write_vec_ave_err(KVP_exc.ave_err());
      KVP_exc_plot.write_polygon([p=pFit,tSep](const double x)->djack_t
				 {
				   return p[0]*exp(-p[2]*(x-tSep))+p[1]*exp(-p[2]*(T-x));
				 }
	,tSep+dt,T-dt,grace::RED);
      // KVP_exc_plot.write_constant_band(0,T,(djack_t)(MP+MP));
      // KVP_exc_plot.write_constant_band(0,T,(djack_t)-MK);
    }
  
  return 0;
}
