#include <glob.h>
#include <tranalisi.hpp>

const size_t T=64,L=32;
const size_t TSep1=10;
const size_t TSep2=12;
const double th=3.9254;
const double pi=th*M_PI/32;

int main()
{
  map<tuple<string,string,string>,vector<double>> data;
  
  char pattern[]="data/*/mes_contr_contr";
  vector<string> conf_list;
  glob_t globbuf;
  if(glob(pattern,0,NULL,&globbuf))
    CRASH("Unable to find pattern %s for conf",pattern);
  
  for(int j=0;j<(int)globbuf.gl_pathc;j++)
    conf_list.push_back(globbuf.gl_pathv[j]);
  
  const size_t nConfs=conf_list.size();
  
  cout<<"found: "<<nConfs<<" confs"<<endl;
  globfree(&globbuf);
  
  /////////////////////////////////////////////////////////////////
  
  for(const auto& conf : conf_list)
    {
      FILE* fin=fopen(conf.c_str(),"r");
      if(fin==nullptr)
	CRASH("%s","");
      
      for(size_t iCombo=0;iCombo<23;iCombo++)
	{
	  char daggered[30],forward[30];
	  const int nProp=fscanf(fin," # Contraction of %s ^ \\dag and %s",daggered,forward);
	  if(nProp==2)
	    for(size_t iBil=0;iBil<5;iBil++)
	      {
		char bil[30];
		const int nBil=fscanf(fin," # %s",bil);
		if(nBil!=1)
		  CRASH("unable to read bil\n");
		
		auto& h=data[{daggered,forward,bil}];
		h.reserve(h.size()+2*T);
		
		for(size_t tri=0;tri<2*T;tri++)
		  {
		    double d;
		    const size_t nTri=fscanf(fin,"%lg\n",&d);
		    if(nTri!=1) CRASH("%zu",nTri);
		    
		    h.push_back(d);
		  }
	      }
	  else
	    CRASH("unable to read props");
	}
      
      fclose(fin);
    }
  
  set_njacks(15);
  
  auto get=
    [&data,&nConfs](const char* back,
		    const char* forw,
		    const char* bil,
		    const size_t ri)
    {
      const auto h=data.find({back,forw,bil});
      if(h==data.end())
	CRASH("Unable to find the correlation %s %s %s",back,forw,bil);
      
      djvec_t res(T);
      
      res=0.0;
      
      jackknivesFill(nConfs,[&res,&h,&ri](const size_t &iConf,
					  const size_t &iClust,
					  const double &w)
      {
	for(size_t t=0;t<T;t++)
	  res[t][iClust]+=w*h->second[ri+2*(t+T*iConf)];
      });
      
      res.clusterize((double)nConfs/njacks);
      
      return res;
    };
  
  auto getP5P5=
    [&get](const char* back,
	   const char* forw)
    {
      return get(back,forw,"P5P5",RE);
    };
  
  // auto getV0P5=
  //   [&get](const char* back,
  // 	   const char* forw)
  //   {
  //     return get(back,forw,"V0P5",RE);
  //   };
  
  auto getVKP5=
    [&get](const char* back,
	   const char* forw)
    {
      return (get(back,forw,"V1P5",IM)+
	      get(back,forw,"V2P5",IM)+
	      get(back,forw,"V3P5",IM))/3;
    };
  
  djack_t Z2_Pi,M_Pi;
  const auto Pi=getP5P5("SMR_L_SMR","SMR_L_SMR");
  two_pts_fit(Z2_Pi,M_Pi,Pi.symmetrized(),T/2,8,T/2,"plots/Pi.xmg");
  
  const auto cPi_ZV_TSep1=get("L_SMR","L_PT1_SMB_L_SMR","V0P5",0);
  const auto cPi_ZV_TSep2=get("L_SMR","L_PT2_SMB_L_SMR","V0P5",0);
  cPi_ZV_TSep1.ave_err().write("plots/corr_Pi_ZV_TSep1.xmg");
  cPi_ZV_TSep2.ave_err().write("plots/corr_Pi_ZV_TSep2.xmg");
  
  cout<<two_pts_corr_fun(Z2_Pi,M_Pi,(size_t)(T/2),T-TSep1,0).ave_err()<<" "<<two_pts_corr_fun(Z2_Pi,M_Pi,(size_t)(T/2),T+TSep1,0)<<endl;
  
  auto fill_Z=
    [](const djvec_t& cZV,
       const djack_t &Z2,
       const djack_t &M,
       const size_t& TSep)
    {
      djvec_t ZV_from_Rat(TSep);
      cout<<" "<<two_pts_corr_fun(Z2,M,(size_t)(T/2),TSep,0).ave_err()<<endl;
      for(size_t t=0;t<TSep;t++) ZV_from_Rat[t]=-0.5*cZV[t]*2*M*exp(M*TSep)/Z2;
      //for(size_t t=TSep+1;t<T;t++) ZV_from_Rat[t]=0.5*cZV[t]/two_pts_corr_fun(Z2,M,(size_t)(T/2),T-TSep,0);
      //ZV_from_Rat[TSep]=ZV_from_Rat[TSep+1]=ZV_from_Rat[T-1]=
      ZV_from_Rat[0]=ZV_from_Rat[TSep/2];
      
      return ZV_from_Rat;
    };
  
  const djvec_t ZV_from_Pi_TSep1=fill_Z(cPi_ZV_TSep1,Z2_Pi,M_Pi,TSep1);
  const djvec_t ZV_from_Pi_TSep2=fill_Z(cPi_ZV_TSep2,Z2_Pi,M_Pi,TSep2);
  ZV_from_Pi_TSep1.ave_err().write("plots/ZV_from_Pi_TSep1.xmg");
  ZV_from_Pi_TSep2.ave_err().write("plots/ZV_from_Pi_TSep2.xmg");
  
  djack_t Z2_B,M_B;
  const auto B=get("SMR_H_SMR","SMR_L_SMR","P5P5",0).symmetrized();
  two_pts_fit(Z2_B,M_B,B,T/2,7,11,"plots/B.xmg");
  cout<<"MB: "<<M_B.ave_err()<<endl;
  
  const auto cB_ZV_TSep1=get("H_SMR","H_PT1_SMB_L_SMR","V0P5",0);
  const auto cB_ZV_TSep2=get("H_SMR","H_PT2_SMB_L_SMR","V0P5",0);
  const djvec_t ZV_from_B_TSep1=fill_Z(cB_ZV_TSep1,Z2_B,M_B,TSep1);
  const djvec_t ZV_from_B_TSep2=fill_Z(cB_ZV_TSep2,Z2_B,M_B,TSep2);
  ZV_from_B_TSep1.ave_err().write("plots/ZV_from_B_TSep1.xmg");
  ZV_from_B_TSep2.ave_err().write("plots/ZV_from_B_TSep2.xmg");
  
  djack_t Z2_D_THP,E_D;
  djack_t Z2_D,M_D;
  const auto D=get("SMC_C_SMC","SMR_L_SMR","P5P5",0).symmetrized();
  const auto D_THP=get("SMC_CTHP_SMC","SMR_L_SMR","P5P5",0).symmetrized();
  two_pts_fit(Z2_D_THP,E_D,D_THP,T/2,7,11,"plots/D_THP.xmg");
  two_pts_fit(Z2_D,M_D,D,T/2,7,11,"plots/D.xmg");
  two_pts_fit(Z2_D_THP,E_D,D_THP,T/2,7,11,"plots/D.xmg");
  cout<<"MD: "<<M_D.ave_err()<<endl;
  const djack_t E_D_cont=cont_en(M_D,pi);
  cout<<"ED: "<<E_D.ave_err()<<endl;
  cout<<"ED_cont: "<<E_D_cont.ave_err()<<endl;
  cout<<"Z2_D: "<<Z2_D.ave_err()<<endl;
  cout<<"Z2_D_THP: "<<Z2_D_THP.ave_err()<<endl;
  
  {
    const double m=0.896;
  
    auto q2=[](const double& mB)
    {
      const double m=0.896;
      const double pi=th*3.14159/32;
      const double e=latt_en(m,pi);
      const double q2=sqr(mB-e)-3*sqr(pi);
      //cout<<"e"<<e<<endl;
      return q2;
    };
    
    cout<<"mB targ: "<<Brent_solve(q2,m,100)<<endl;
  }
  
  const djack_t Q2max=sqr(M_B-M_D);
  const djack_t Q2=sqr(M_B-E_D)-3*sqr(th*M_PI/L);
  const djack_t Q2_cont=sqr(M_B-E_D_cont)-3*sqr(th*M_PI/L);
  cout<<"Q2max: "<<Q2max.ave_err()<<endl;
  cout<<"Q2: "<<Q2.ave_err()<<endl;
  cout<<"Q2_cont: "<<Q2_cont.ave_err()<<endl;
  
  const djvec_t cdE=effective_mass(B)-effective_mass(D_THP);
  const djack_t dE=constant_fit(cdE,4,8,"plots/dE.xmg");
  
  const djvec_t cBtoD_VK_TSep1=(getVKP5("H_SMR","CTHP_PT1_SMD_L_SMR")-getVKP5("H_SMR","CTHM_PT1_SMD_L_SMR"))/2;
  const djvec_t cBtoD_VK_TSep2=(getVKP5("H_SMR","CTHP_PT2_SMD_L_SMR")-getVKP5("H_SMR","CTHM_PT2_SMD_L_SMR"))/2;
  // effective_mass(cBtoD_V0_TSep1,T/2,0).ave_err().write("plots/BtoD_V0_TSep1.xmg");
  // effective_mass(cBtoD_V0_TSep2,T/2,0).ave_err().write("plots/BtoD_V0_TSep2.xmg");
  djvec_t BtoD_VK_TSep1(TSep1);
  djvec_t BtoD_VK_TSep2(TSep2);
  for(size_t t=0;t<TSep1;t++)
    BtoD_VK_TSep1[t]=-cBtoD_VK_TSep1[t]*exp(M_B*t+E_D*(TSep1-t));
  for(size_t t=0;t<TSep2;t++)
    BtoD_VK_TSep2[t]=-cBtoD_VK_TSep2[t]*exp(M_B*t+E_D*(TSep2-t));
  
  BtoD_VK_TSep1.ave_err().write("plots/BtoD_VK_TSep1.xmg");
  BtoD_VK_TSep2.ave_err().write("plots/BtoD_VK_TSep2.xmg");
  
  grace_file_t BtoD_V0_plot("plots/BtoD_V0.xmg");
  BtoD_V0_plot.set_no_line();
  BtoD_V0_plot.set_all_colors(grace::BLACK);
  BtoD_V0_plot.set_legend("tsep=10");
  for(size_t t=0;t<TSep1;t++)
    BtoD_V0_plot.write_ave_err(t/(double)TSep1,BtoD_VK_TSep1[t].ave_err());
  
  BtoD_V0_plot.new_data_set();
  BtoD_V0_plot.set_no_line();
  for(size_t t=0;t<TSep2;t++)
    BtoD_V0_plot.write_ave_err(t/(double)TSep2,BtoD_VK_TSep2[t].ave_err());
  BtoD_V0_plot.set_legend("tsep=12");
  BtoD_V0_plot.set_all_colors(grace::RED);
  BtoD_V0_plot.close();
  
  const auto cDth_ZV_TSep1=get("CTHM_SMC","CTHP_PT1_SMD_L_SMR","V0P5",0);
  const auto cDth_ZV_TSep2=get("CTHP_SMC","CTHM_PT2_SMD_L_SMR","V0P5",0);
  const auto cDthP_VK_TSep2=getVKP5("CTHP_SMC","CTHM_PT2_SMD_L_SMR");
  const auto cDthM_VK_TSep2=getVKP5("CTHM_SMC","CTHP_PT2_SMD_L_SMR");
  const djvec_t ZV_from_Dth_TSep1=fill_Z(cDth_ZV_TSep1,Z2_D_THP,E_D,TSep1);
  const djvec_t ZV_from_Dth_TSep2=fill_Z(cDth_ZV_TSep2,Z2_D_THP,E_D,TSep2);
  ZV_from_Dth_TSep1.ave_err().write(vector_up_to<double>(TSep1)/TSep1,"plots/ZV_from_Dth_TSep1.xmg");
  ZV_from_Dth_TSep2.ave_err().write(vector_up_to<double>(TSep2)/TSep2,"plots/ZV_from_Dth_TSep2.xmg");
  cDthP_VK_TSep2.ave_err().write(vector_up_to<double>(TSep2)/TSep2,"plots/cDthP_VK_TSep2.xmg");
  cDthM_VK_TSep2.ave_err().write("plots/cDthM_VK_TSep2.xmg");
  effective_mass(cDthM_VK_TSep2,T/2,0).ave_err().write("/tmp/o.xmg");
  const auto cD_ZV_TSep1=get("C_SMC","C_PT1_SMD_L_SMR","V0P5",0);
  //const auto cD_ZV_TSep2=get("C_SMC","C_PT2_SMD_L_SMR","V0P5",0);
  const djvec_t ZV_from_D_TSep1=fill_Z(cD_ZV_TSep1,Z2_D,M_D,TSep1);
  //const djvec_t ZV_from_D_TSep2=fill_Z(cD_ZV_TSep2,Z2_D,M_D,TSep2);
  ZV_from_D_TSep1.ave_err().write(vector_up_to<double>(TSep1)/TSep1,"plots/ZV_from_D_TSep1.xmg");
  //ZV_from_D_TSep2.ave_err().write("plots/ZV_from_D_TSep2.xmg");
  
  return 0;
}
