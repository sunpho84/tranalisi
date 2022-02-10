#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "effective.hpp"
#include "fit.hpp"
#include "gevp.hpp"
#include "meas_vec.hpp"

#include <GM3/kernels.hpp>
#include <GM3/params.hpp>
#include <GM3/perens.hpp>
#include <GM3/VKVKRepresentation.hpp>

void perens_t::computeAmu(const RegoType& rego,
			  const jack_t<VKVKRepFiniteVol>& rep,
			  const jack_t<VKVKRepInfiniteVol>& infVolRep)
{
  using namespace Eigen;
  
    // for(const RegoType& rego : {REGO_TM,REGO_OS})
  //   {
  console<<"Computing amu fore rego "<<regoTag[rego]<<endl;
  
  const djack_t aMRho=
    rep.actOn([](const VKVKRepFiniteVol& rep)
    {
     return
       rep.mRho;
    });
  
  const djack_t aFromMRho=
    aMRho/0.775;
  
  const djack_t Zrego=
    Z[regoZId[rego]];
  
  // vector<djvec_t> eig;
  // vector<djvec_t> recastEigvec;
  // vector<djvec_t> origEigvec;
  
  // const djvec_t c00=getAveForRego(0,nSources,1,rego);
  // const djvec_t c01=-getAveForRego(0,nSources,4,rego);
  // const djvec_t c11=-getAveForRego(0,nSources,2,rego);
  // const size_t t0=3;
  
  // tie(eig,recastEigvec,origEigvec)=gevp({c00,c01,c01,c11},t0);
  
  // eig[0].ave_err().write("plots/eig1Rego"+regoTag[rego]+".xmg");
  // eig[1].ave_err().write("plots/eig2Reno"+regoTag[rego]+".xmg");
  
  // const djack_t eig0MDiagFit=constant_fit(effective_mass(eig[0]),tMinFit,tMaxFit,"plots/eff_eig1.xmg");
  // console<<"eig0 mass: "<<eig0MDiagFit.ave_err()<<endl;
  
  // const djack_t eig1MDiagFit=constant_fit(effective_mass(eig[1]),13,18,"plots/eff_eig2.xmg");
  // out<<"eig1 mass: "<<eig1MDiagFit.ave_err()<<endl;
  
  // djvec_t SL0(THp1),SS0(THp1);
  // djvec_t SL1(THp1),SS1(THp1);
  // for(size_t t=0;t<=TH;t++)
  // 	{
  // 	  for(size_t ijack=0;ijack<=njacks;ijack++)
  // 	    {
  // 	      typedef Matrix<double,2,2> Matr;
  
  // 	      Matr e;
  // 	      const auto& ei=origEigvec;
  // 	      e(0,0)=ei[0][t][ijack];
  // 	      e(0,1)=ei[2][t][ijack];
  // 	      e(1,0)=ei[1][t][ijack];
  // 	      e(1,1)=ei[3][t][ijack];
  
  // 	      Matr c;
  // 	      c(0,0)=c00[t][ijack];
  // 	      c(0,1)=c01[t][ijack];
  // 	      c(1,0)=c01[t][ijack];
  // 	      c(1,1)=c11[t][ijack];
  
  // 	      const Matr sl=c*e.transpose();
  // 	      const Matr ss=e*c*e.transpose();
  
  // 	      SL0[t][ijack]=sl(0,0);
  // 	      SS0[t][ijack]=ss(0,0);
  // 	      SL1[t][ijack]=sl(0,1);
  // 	      SS1[t][ijack]=ss(1,1);
  // 	    }
  // 	}
  
  // djack_t eig0ZS,eig0ZL,eig0M;
  // two_pts_SL_fit(eig0ZS,eig0ZL,eig0M,SL0,SS0,TH,tMinVKVK,tMaxVKVK,"plots/SL0Rego"+regoTag[rego]+".xmg");
  // const djack_t eig0Z2L=eig0ZL*eig0ZL;
  // console<<"Z20L: "<<eig0Z2L<<endl;
  
  // djvec_t subCorr1=c00;
  // for(size_t t=0;t<=TH;t++)
  // 	subCorr1[t]-=two_pts_corr_fun(eig0Z2L,eig0M,TH,t,+1);
  
  // djack_t eig1Z2L,eig1M;
  // two_pts_fit(eig1Z2L,eig1M,subCorr1,TH,8,19,"plots/SL1Rego"+regoTag[rego]+".xmg");
  // console<<"Z21L: "<<eig1Z2L<<endl;
  // console<<"M1L: "<<eig1M<<endl;
  //getAveForRego(0,nSources,1,rego).bin_write("plots/VKVK_corr_fromAmu_"+regoTag[rego]+".dat");
  getAveForRego(0,nSources,1,rego).ave_err().write("plots/VKVK_corr_fromAmu_"+regoTag[rego]+".xmg");
  
  const double eu=2.0/3,ed=-1.0/3;
  const djvec_t corr=
    getAveForRego(0,nSources,1,rego)*sqr(Zrego)*(sqr(eu)+sqr(ed));
  
  djack_t mVK1,Z2VK1;
  two_pts_fit(Z2VK1,mVK1,corr,TH,tMinVKVK,tMaxVKVK,"plots/eff_mass_VKVK_twopts_fit_rego"+regoTag[rego]+".xmg");
  // console<<"Z2: "<<Z2VK1<<endl;
  // // djack_t mVK2,Z2VK2;
  // // two_pts_fit(Z2VK2,mVK2,corr,TH,22,32);
  
  grace_file_t amu("plots/amuRego"+regoTag[rego]+".xmg");
  djvec_t amuInt(TH),amuSubs1(TH),amuSubs2(TH),amuSubs3(TH),amuSubs4(TH),amuSubs5(TH),amuSubs6(TH);
  for(size_t upto=0;upto<TH;upto++)
    {
      djvec_t corrRefattaSingleExp=corr;
      djvec_t corrRefattaTFinito=corr;
      djvec_t corrRefattaTInfinito=corr;
      for(size_t t=upto;t<=TH;t++)
	{
	  corrRefattaSingleExp[t]=two_pts_corr_fun(Z2VK1,mVK1,TH,t,0);
	  corrRefattaTFinito[t]=rep(t)();
	  corrRefattaTInfinito[t]=infVolRep(t)();
	}
      
      size_t THm1=TH-1;
      const djack_t ja=a;
      const djack_t jal=aFromPhysLine;
      const djack_t jam=aFromMRho;
      const djack_t cInt=integrate_corr_times_kern_up_to(corr,T,ja,upto)*1e10;
      const djack_t cSubs1=integrate_corr_times_kern_up_to(corrRefattaSingleExp,T,ja,THm1)*1e10;
      const djack_t cSubs2=integrate_corr_times_kern_up_to(corrRefattaTFinito,T,ja,THm1)*1e10;
      const djack_t cSubs3=integrate_corr_times_kern_up_to(corrRefattaTFinito,T,jal,THm1)*1e10;
      const djack_t cSubs4=integrate_corr_times_kern_up_to(corrRefattaTInfinito,T,jal,THm1)*1e10;
      const djack_t cSubs5=integrate_corr_times_kern_up_to(corrRefattaTFinito,T,jam,THm1)*1e10;
      const djack_t cSubs6=integrate_corr_times_kern_up_to(corrRefattaTInfinito,T,jam,THm1)*1e10;
      amuInt[upto]=cInt;
      amuSubs1[upto]=cSubs1;
      amuSubs2[upto]=cSubs2;
      amuSubs3[upto]=cSubs3;
      amuSubs4[upto]=cSubs4;
      amuSubs5[upto]=cSubs5;
      amuSubs6[upto]=cSubs6;
      
      if(upto==0 and rego==REGO_TM)
	{
	  ofstream outForNazOri("/tmp/B07264.dat");
	  outForNazOri.precision(16);
	  outForNazOri<<njacks<<" "<<T<<endl;
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    for(size_t t=0;t<T;t++)
	      outForNazOri<<t<<" "<<((t==0)?0:corr[(t<TH)?t:(T-t)][ijack])/(sqr(Zrego[ijack])*(sqr(eu)+sqr(ed)))<<endl;
	  
	  ofstream outForNazFilt("/tmp/B07264_filtered.dat");
	  outForNazFilt.precision(16);
	  outForNazFilt<<njacks<<" "<<T<<endl;
	  for(size_t ijack=0;ijack<njacks;ijack++)
	    for(size_t t=0;t<T;t++)
	      outForNazFilt<<t<<" "<<((t==0)?0:corrRefattaTFinito[(t<TH)?t:(T-t)][ijack])/(sqr(Zrego[ijack])*(sqr(eu)+sqr(ed)))<<endl;
	}
    }
  amu.set_xaxis_label("t");
  
  djvec_t corrRefatta(THp1);
  for(size_t t=0;t<TH;t++)
    corrRefatta[t]=
      rep(t)();// two_pts_corr_fun(eig0Z2L,eig0M,TH,t,0)+
  // two_pts_corr_fun(eig1Z2L,eig1M,TH,t,0);
  corrRefatta.ave_err().write("plots/corr_VKVK_refatta"+regoTag[rego]+".xmg");
  corrRefatta.bin_write("plots/corr_VKVK_refatta"+regoTag[rego]+".dat");
  effective_mass(corrRefatta,TH,0).ave_err().write("plots/eff_mass_VKVK_refatta_Rego"+regoTag[rego]+".xmg");
  
  amu.write_vec_ave_err(amuInt.ave_err());
  amu.set_no_line();
  amu.set_legend("Pure integration");
  amu.set_all_colors(grace::RED);
  
  amu.write_vec_ave_err(amuSubs1.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::BLUE);
  amu.set_legend("Analytic continuation of ground state");
  
  amu.write_vec_ave_err(amuSubs2.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::ORANGE);
  amu.set_legend("A la luscher rep");
  
  amu.write_vec_ave_err(amuSubs3.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::ORANGE);
  amu.set_legend("A la luscher rep + lt");
  
  amu.write_vec_ave_err(amuSubs4.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::MAGENTA);
  amu.set_legend("A la luscher rep + lt + infvol");
  
  amu.write_vec_ave_err(amuSubs5.ave_err());
  amu.set_no_line();
  amu.set_all_colors(grace::MAROON);
  amu.set_legend("A la luscher rep + lt rho + infvol");
  
  amu.new_data_set();
  amu.set_legend("BMW light connected");
  amu.set_all_colors(grace::GREEN4);
  amu.write_constant_band(0,TH,djack_t(gauss_filler_t{652.4,6.0,23423}));
  
  console<<"amu: "<<(aAve*aAve)<<" "<<amuSubs3[20]<<endl;
  console<<"amu: "<<(aAve*aAve)<<" "<<amuSubs4[6]<<" inf vol"<<endl;
  console<<"amu from rho: "<<(aAve*aAve)<<" "<<amuSubs5[6]<<endl;
  console<<"amu from rho: "<<(aAve*aAve)<<" "<<amuSubs6[6]<<" inf vol"<<endl;
  
  const djack_t ddss=sqrt(a/aFromMRho);
  console<<"amu from rho with prefact: "<<(aAve*aAve)<<" "<<amuSubs5[6]*ddss<<endl;
}

