#include "effective.hpp"
#include "fit.hpp"
#include "grace.hpp"
#include "meas_vec.hpp"
#include "raw_file.hpp"
#include "tools.hpp"
#include <functional>
#include <tranalisi.hpp>
#include <vector>

size_t L;
size_t T;
double TH;
double zp_fr_zs;
double aml;
size_t tmin;
size_t tmax;

template <typename T,
	  typename TV>
T fitSlope2Ansatz(const TV& pars,const T&aM,double t)
{
  T th=tanh(aM*(t-TH));
  return pars[0]+pars[1]*sqr(t-TH)+pars[2]*(t-TH)*th;
}

vector<complex<double>> load_single_bubble(const string& TAG)
{
  cout<<"======== loading bubble "<<TAG<<" =========="<<endl;
  
  vector<complex<double>> d;
  
  ifstream is("jacks/bubble_"+TAG+"_ingr");
  
  double r,i;
  while(is>>r>>i)
    d.push_back({r,i});
  
  return d;
}

djvec_t load_bubble(const string& TAG1,const string& TAG2,const bool useRe,const bool useIm,const bool sub=true)
{
  vector<complex<double>> d1=load_single_bubble(TAG1);
  
  vector<complex<double>> d2=(TAG1!=TAG2)?load_single_bubble(TAG2):d1;
  
  if(d1.size()!=d2.size())
    CRASH("Number of confs not matching, %zu vs %zu",d1.size(),d2.size());
  
  const int nconfs_max=d1.size()/T;
  cout<<"NConfs: "<<nconfs_max<<endl;
  
  const int clust_size=nconfs_max/njacks;
  const int nconfs=clust_size*njacks;
  if(nconfs!=nconfs_max)
    cout<<"NConfs adjusted: "<<nconfs<<endl;
  
  // complex<double> ave={0.0,0.0};
  // for(int iconf=0;iconf<nconfs;iconf++)
  //   for(int dt=0;dt<(int)T;dt++)
  //     ave+=d[dt+T*iconf];
  // ave/=T*nconfs;
  // for(int iconf=0;iconf<nconfs;iconf++)
  //   for(int dt=0;dt<(int)T;dt++)
  //     d[dt+T*iconf]-=ave;
  
  djvec_t aveBubble1(2),aveBubble2(2);
  djack_t aveBubbleSquare;
  djvec_t avePerTimeRe1(T),avePerTimeRe2(T);
  djvec_t avePerTimeIm1(T),avePerTimeIm2(T);
   for(int iconf=0;iconf<nconfs;iconf++)
     {
       const int iclust=iconf/clust_size;
      
       for(int dt=0;dt<(int)T;dt++)
	 {
	  avePerTimeRe1[dt][iclust]+=d1[dt+T*iconf].real();
	  avePerTimeIm1[dt][iclust]+=d1[dt+T*iconf].imag();
	  avePerTimeRe2[dt][iclust]+=d2[dt+T*iconf].real();
	  avePerTimeIm2[dt][iclust]+=d2[dt+T*iconf].imag();
	   
	  double r1=d1[dt+T*iconf].real();
	  double i1=d1[dt+T*iconf].imag();
	  double r2=d2[dt+T*iconf].real();
	  double i2=d2[dt+T*iconf].imag();
	   
	   if(useRe)
	     {
	       aveBubbleSquare[iclust]+=r1*r2;
	       aveBubble1[RE][iclust]+=r1;
	       aveBubble2[RE][iclust]+=r2;
	     }
	   
	   if(useIm)
	     {
	       aveBubbleSquare[iclust]+=i1*i2;
	       aveBubble1[IM][iclust]+=i1;
	       aveBubble2[IM][iclust]+=i2;
	     }
	 }
     }
   aveBubble1.clusterize(clust_size);
   aveBubble2.clusterize(clust_size);
   aveBubbleSquare.clusterize(clust_size);
   djack_t aveBubbleSub;
   if(useRe)
     aveBubbleSub+=aveBubble1[RE]*aveBubble2[RE];
   if(useIm)
     aveBubbleSub+=aveBubble1[IM]*aveBubble2[IM];
   
   aveBubbleSub-=aveBubbleSquare;
   aveBubbleSub/=(T*(T-1));
   
   avePerTimeRe1.clusterize(clust_size);
   avePerTimeIm1.clusterize(clust_size);
   avePerTimeRe2.clusterize(clust_size);
   avePerTimeIm2.clusterize(clust_size);
   
   // aveBubble/=T;
   // aveBubbleSquare/=T;
   
   // double test=0;
   // int n=0;
   // for(int iconf=0;iconf<nconfs;iconf++)
   //   for(int jconf=0;jconf<nconfs;jconf++)
   //     if(iconf!=jconf)
   // 	 for(int t1=0;t1<(int)T;t1++)
   // 	   for(int t2=0;t2<(int)T;t2++)
   // 	     if(t1!=t2)
   // 	       {
   // 		 test+=
   // 		   (d[t1+T*iconf]*conj(d[t2+T*jconf])).real();
   // 		 n++;
   // 	       }
  // aveBubble.clusterize(clust_size);
  // aveBubble/=T;
  // cout<<"Ave of bubble: "<<aveBubble.ave_err()<<endl;
  
  djvec_t bubbleProd(T);
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      const int iclust=iconf/clust_size;
      
      vec_ave_err_t single_conf_bubble(T);
      
      for(int dt=0;dt<(int)T;dt++)
	{
	  double sx=0.0,s2x=0.0;
	  
	  for(int t1=0;t1<(int)T;t1++)
	    {
	      int t2=(t1+dt)%T;
	      double x=0.0;
	      if(useRe)
		x=d1[t1+T*iconf].real()*d2[t2+T*iconf].real();
	      if(useIm)
		x+=d1[t1+T*iconf].imag()*d2[t2+T*iconf].imag();
	      bubbleProd[dt][iclust]+=x;
	      
	      sx+=x;
	      s2x+=x*x;
	    }
	  
	  sx/=T;
	  s2x/=T;
	  s2x-=sx*sx;
	  
	  single_conf_bubble[dt]={sx,sqrt(s2x/(T-1))};
	}
      
       if(iconf<=10)
	single_conf_bubble.write("plots/singleConf"+TAG1+to_string(iconf)+".xmg");
    }
  
  bubbleProd.clusterize(clust_size);
  bubbleProd/=T;
  
  djvec_t aveBubbleProdTime(T);
  for(int dt=0;dt<(int)T;dt++)
    for(int t1=0;t1<(int)T;t1++)
      {
	int t2=(t1+dt)%T;
	if(RE)
	  aveBubbleProdTime[dt]+=avePerTimeRe1[t1]*avePerTimeRe2[t2];
	if(IM)
	  aveBubbleProdTime[dt]+=avePerTimeIm1[t1]*avePerTimeIm2[t2];
      }
  aveBubbleProdTime/=T;
  //const djack_t aveBubbleProd=(aveBubbleProdTime.sum()-aveBubbleProdTime[0])/(T-1);
  ((djvec_t)((djvec_t)(bubbleProd-aveBubbleProdTime)).symmetrized()*pow(L,3)).ave_err().write("plots/aveBubbleProd.xmg");
  
  cout.precision(16);
  cout<<"Bubble[1]: "<<bubbleProd[1].ave_err()<<endl;
  cout<<"Bubble[2]: "<<bubbleProd[2].ave_err()<<endl;
  cout<<"Bubble["<<TH-1<<"]: "<<bubbleProd[TH-1].ave_err()<<endl;
  const djack_t bias=aveBubble1[IM]*aveBubble2[IM]/sqr(T)-aveBubbleSub;
  cout<<"Ave: "<<aveBubbleSub.ave_err()<<" biassed: "<<((djack_t)(aveBubble1[IM]*aveBubble2[IM]/sqr(T))).ave_err()<<" bias: "<<bias.ave_err()<<endl;
  //cout<<"Time averaged: "<<aveBubbleProd.ave_err()<<endl;
  // cout<<"Test: "<<test/n<<endl;
  
  const djack_t y=bubbleProd[2]/aveBubbleSub-1;
  cout<<"Rel diff: "<<y.ave_err()<<endl;
  
  if(sub)
    bubbleProd-=
      aveBubbleSub;
   
   // bubbleProd-=(djack_t)bubbleProd[TH];
   
   
  
  // bubbleProd-=aveBubble;
  
  // vector<complex<double>> dPerConf(nconfs);
  // for(int iconf=0;iconf<nconfs;iconf++)
  //   {
  //     dPerConf[iconf]=0;
  //     for(int t=0;t<(int)T;t++)
  // 	dPerConf[iconf]+=d[t+T*iconf];
  //     dPerConf[iconf]/=T;
  //   }
  
  // djack_t a,a2;
  // for(int ijack=0;ijack<=(int)njacks;ijack++)
  //   {
  //     int n=0;
  //     a[ijack]=0;
  //     int iMinJConf=ijack*clust_size;
  //     int iMaxJConf=iMinJConf+clust_size;
  //     for(int iconf=0;iconf<nconfs;iconf++)
  // 	if(iconf<iMinJConf or iconf>iMaxJConf)
  // 	  for(int jconf=0;jconf<nconfs;jconf++)
  // 	    if(jconf!=iconf and (jconf<iMinJConf or jconf>iMaxJConf))
  // 	      {
  // 		// if(RI==0)
  // 		  a[ijack]+=dPerConf[iconf].real()*dPerConf[jconf].real();
  // 		// else
  // 		  a[ijack]+=dPerConf[iconf].imag()*dPerConf[jconf].imag();
  // 		n++;
  // 	      }
      
  //     a[ijack]/=n;
  //   }
  
  //bubble-=a;
  
  // cout<<"Other average: "<<a.ave_err()<<endl;
  
  return bubbleProd.symmetrized()*pow(L,3);
}

// djvec_t load_bubble_P5P5()
// {
//   ifstream is("jacks/bubble_P5_B1_ingr");
  
//   vector<double> d;
//   double r,i;
//   while(is>>r>>i)
//     d.push_back(i);
  
//   const int nconfs_max=d.size()/T;
//   cout<<"NConfs: "<<nconfs_max<<endl;
  
//   const int clust_size=nconfs_max/njacks;
//   const int nconfs=clust_size*njacks;
//   if(nconfs!=nconfs_max)
//     cout<<"NConfs adjusted: "<<nconfs<<endl;
  
//   djvec_t avePerTime(T);
//   for(int iconf=0;iconf<nconfs;iconf++)
//     {
//       const int iclust=iconf/clust_size;
      
//       for(int dt=0;dt<(int)T;dt++)
// 	avePerTime[dt][iclust]+=d[dt+T*iconf];
//     }
//   avePerTime.clusterize(clust_size);
  
//   djvec_t bubble(T);
//   for(int iconf=0;iconf<nconfs;iconf++)
//     {
//       const int iclust=iconf/clust_size;
      
//       for(int dt=0;dt<(int)T;dt++)
// 	for(int t1=0;t1<(int)T;t1++)
// 	  {
// 	    int t2=(t1+dt)%T;
// 	    bubble[dt][iclust]+=d[t1+T*iconf]*d[t2+T*iconf];
// 	  }
//     }
  
//   bubble.clusterize(clust_size);
//   bubble/=T;
  
//   djvec_t aveBubble(T);
//   for(int dt=0;dt<(int)T;dt++)
//     for(int t1=0;t1<(int)T;t1++)
//       {
// 	int t2=(t1+dt)%T;
// 	aveBubble[dt]+=avePerTime[t1]*avePerTime[t2];
//       }
//   aveBubble/=T;
//   //cout<<aveBubble.ave_err()<<endl;
//   bubble-=aveBubble;
  
//   return bubble.symmetrized()*pow(L,3);
// }


int main()
{
  input_file_t input("pars.txt");
  set_njacks(input.read<size_t>("njacks"));
  L=input.read<size_t>("L");
  T=input.read<size_t>("T");
  tmin=input.read<size_t>("tmin");
  tmax=input.read<size_t>("tmax");
  TH=T/2.0;
  aml=input.read<double>("aml");
  zp_fr_zs=input.read<double>("zp_fr_zs");
  
  const djvec_t P5P5_TM_SL=read_djvec("jacks/P5P5_TM",T,0).symmetrized();
  P5P5_TM_SL.ave_err().write("plots/corr_P5P5_TM_SL.xmg");
  const djvec_t P5P5_TM_SS=read_djvec("jacks/P5P5_TM_sme",T,0).symmetrized();
  const djvec_t eff_P5P5_TM_SS=effective_mass(P5P5_TM_SS);
  P5P5_TM_SS.ave_err().write("plots/corr_P5P5_TM_SS.xmg");
  
  const djvec_t P5P5_OS_SL=read_djvec("jacks/P5P5_OS",T,0).symmetrized();
  P5P5_OS_SL.ave_err().write("plots/corr_P5P5_OS_SL.xmg");
  const djvec_t P5P5_OS_SS=read_djvec("jacks/P5P5_OS_sme",T,0).symmetrized();
  P5P5_OS_SS.ave_err().write("plots/corr_P5P5_TM_SS.xmg");
  
  djack_t aM_TM,ZS_TM,ZL_TM;
  two_pts_SL_fit(ZS_TM,ZL_TM,aM_TM,P5P5_TM_SL,P5P5_TM_SS,TH,tmin,tmax,"plots/eff_mass_P5P5_TM.xmg");
  djack_t aM_OS,ZS_OS,ZL_OS;
  two_pts_SL_fit(ZS_OS,ZL_OS,aM_OS,P5P5_OS_SL,P5P5_OS_SS,TH,tmin,tmax,"plots/eff_mass_P5P5_OS.xmg");
  const djack_t af_TM=2*aml*ZL_TM/sqr(aM_TM);
  
  cout<<"aM_TM: "<<aM_TM.ave_err()<<endl;
  cout<<"af_TM: "<<af_TM.ave_err()<<endl;
  
  cout<<"aM_OS: "<<aM_OS.ave_err()<<endl;
  
  
  const djack_t ZP_fr_ZS=ZL_OS/ZL_TM;
  cout<<"Zp/Zs: "<<ZP_fr_ZS.ave_err()<<endl;
  
  zp_fr_zs=ZP_fr_ZS.ave();
  
  const djack_t norm=-sqr(af_TM)*sqr(aml)/pow(aM_TM,3)/sqr(zp_fr_zs);
  cout<<"Norm: "<<norm.ave_err()<<endl;
  
  // const djvec_t P5P5_disco=load_bubble_P5P5();
  // P5P5_disco.ave_err().write("plots/P5P5_disco.xmg");
  // const djvec_t eff_P5P5_disco=effective_mass(P5P5_disco);
  // eff_P5P5_disco.ave_err().write("plots/eff_P5P5_disco.xmg");
  
  
  const djvec_t P5P5_SS=read_djvec("jacks/P5P5_SS_TM",T,0).symmetrized()*norm;
  P5P5_SS.ave_err().write("plots/P5P5_SS_TM.xmg");
  
  const djvec_t rat_conn=P5P5_SS/P5P5_TM_SS;
  rat_conn.ave_err().write("plots/ratio_conn_TM.xmg");
  const djvec_t eff_slope_conn=effective_slope(rat_conn,eff_P5P5_TM_SS,TH);
  eff_slope_conn.ave_err().write("plots/eff_slope_conn_TM.xmg");
  
  // djvec_t rat_conn_sub=rat_conn;
  // for(size_t t=0;t<=TH;t++)
  //   {
  //     double dt=(double)TH-t;
  //     rat_conn_sub[t]-=sqr(SL)*sqr(dt);
  //   }
  // rat_conn_sub.ave_err().write("plots/ratio_conn_sub_TM.xmg");
  // const djvec_t eff_slope_conn_sub=effective_slope(rat_conn_sub,eff_P5P5,TH);
  // eff_slope_conn_sub.ave_err().write("plots/eff_slope_conn_sub_TM.xmg");
  // const djvec_t eff_slope_conn_sub_sub=eff_slope_conn_sub-2*A1_fr_A*SL;
  // const djack_t conn_l7_fit=constant_fit(eff_slope_conn_sub_sub,tmin,TH-1,"plots/eff_slope_conn_sub_sub_TM.xmg");
  // const djack_t l7=conn_l7_fit*norm;
  // cout<<"l7: "<<l7.ave_err()<<endl;
  
  const djvec_t S0P5_S0P5=load_bubble("S0P5_TM","S0P5_TM",true,true,false)*norm;
  S0P5_S0P5.ave_err().write("plots/S0P5_S0P5_TM.xmg");
  
  const djvec_t rat_disco=S0P5_S0P5/P5P5_TM_SS;
  rat_disco.ave_err().write("plots/ratio_disco.xmg");
  const djvec_t eff_slope_disco=effective_slope(rat_disco,eff_P5P5_TM_SS,TH);
  eff_slope_disco.ave_err().write("plots/eff_slope_disco_TM.xmg");
  
  const djvec_t rat_combo=rat_conn-rat_disco;
  const djvec_t eff_slope_combo=effective_slope(rat_combo,eff_P5P5_TM_SS,TH);
  const djack_t l7=constant_fit(eff_slope_combo,5,10,"plots/eff_slope_combo.xmg");
  cout<<"l7: "<<l7.ave_err()<<endl;
  
  auto fit=[](const djack_t& aM,const djvec_t& eff_P5P5,const djvec_t rat,const string tag="")
  {
    cout<<"FIT OF: "<<tag<<endl;
    cout<<"===================="<<endl;
    
    jack_fit_t fitter;
    djvec_t pars(3);
    fitter.add_fit_par(pars[0],"C",1000,0.1);
    fitter.add_fit_par(pars[1],"Sl",0,0.1);
    fitter.add_fit_par(pars[2],"E",0,0.1);
    
    const size_t tmax=TH;
    for(size_t it=tmin;it<=tmax;it++)
      fitter.add_point(rat[it],[t=(double)it,aM](const vector<double>& pars,const size_t iel)
      {
	return fitSlope2Ansatz(pars,aM[iel],t);
      });
    fitter.fit();
    
    cout<<pars.ave_err()<<endl;
    
    auto effSlope=[aM,pars](const double t)->djack_t
    {
      djvec_t pars0=pars;
      pars0[1]=0.0;
      pars0[2]=1.0;
      
      return
	(fitSlope2Ansatz(pars,aM,t+1)-fitSlope2Ansatz(pars,aM,t))/
	(fitSlope2Ansatz(pars0,aM,t)-fitSlope2Ansatz(pars0,aM,t+1));
    };
    
    grace_file_t ratio_fit("plots/fit_ratio_"+tag+".xmg");
    ratio_fit.write_vec_ave_err(rat.ave_err());
    ratio_fit.write_polygon([&pars,&aM](const double t){return fitSlope2Ansatz(pars,aM,t);},tmin,tmax);
    
    grace_file_t eff_slope_fit("plots/eff_slope_"+tag+".xmg");
    eff_slope_fit.write_vec_ave_err(effective_slope(rat,eff_P5P5,TH).ave_err());
    eff_slope_fit.write_polygon(effSlope,tmin,tmax);
    
    return pars;
  };
  
  const djvec_t pars_conn=fit(aM_TM,eff_P5P5_TM_SS,rat_conn,"conn");
  const djvec_t pars_disco=fit(aM_TM,eff_P5P5_TM_SS,rat_disco,"disco");
  
  // const djvec_t pars_combo=fit(aM,eff_P5P5,rat_combo*norm,"combo");
  // const djack_t l7combo=pars_combo[2];
  // const djack_t l7sub=l7combo-T*pars_combo[1];
  // cout<<"l7 combo: "<<l7combo.ave_err()<<endl;
  // cout<<"l7 sub: "<<l7sub.ave_err()<<endl;
  
  // const djack_t estim_zp_fr_zs=sqrt(pars_disco[1]/pars_conn[1]);
  // const djack_t fact=estim_zp_fr_zs/zp_fr_zs;
  // cout<<"Estimator: "<<fact.ave_err()<<endl;
  
  // const djvec_t rat_eff_combo=norm*(rat_conn-rat_disco/sqr(estim_zp_fr_zs));
  // const djvec_t pars_eff_combo=fit(aM,eff_P5P5,rat_eff_combo,"eff_combo");
  // const djack_t l7_eff_combo=pars_eff_combo[2];
  // cout<<"l7 eff combo: "<<l7_eff_combo.ave_err()<<endl;
  
  // const djack_t l7_eff_combo_sub=l7_eff_combo-pars_eff_combo[1]*T;
  // cout<<"l7 eff sub: "<<l7_eff_combo_sub.ave_err()<<endl;
  
      // const djvec_t P5_B1=load_bubble("P5_B1",false,true);
      // P5_B1.ave_err().write("plots/P5_B1.xmg");
      // effective_mass(P5_B1).ave_err().write("plots/eff_mass_P5_B1.xmg");
      
      // const djvec_t P5_B1_eta=load_bubble("P5_B1",true,false,false);
      // P5_B1_eta.ave_err().write("plots/P5_B1_eta.xmg");
      // effective_mass(P5_B1_eta).ave_err().write("plots/eff_mass_P5_B1_eta.xmg");
      
      // const djvec_t P5P5_OS=read_djvec("jacks/P5P5_OS",T,0).symmetrized();
      // P5P5_OS.ave_err().write("plots/P5P5_OS_corr.xmg");
      // const djvec_t eff_P5P5_OS=effective_mass(P5P5_OS);
      // const djack_t aM_OS=constant_fit(eff_P5P5_OS,tmin,TH,"plots/eff_mass_P5P5_OS.xmg");
      
      // const djvec_t P5P5_pi0=2*P5_B1+P5P5_OS;
      // P5P5_pi0.ave_err().write("plots/P5P5_pi0_corr.xmg");
      // const djvec_t eff_P5P5_pi0=effective_mass(P5P5_pi0);
      // const djack_t aM_pi0=constant_fit(eff_P5P5_pi0,tmin,TH,"plots/eff_mass_P5P5_pi0.xmg");
      
      // const djvec_t P5P5_pi0_sub=P5P5_pi0-P5P5_pi0[TH]+P5P5_pi0[TH].ave();
      // P5P5_pi0_sub.ave_err().write("plots/P5P5_pi0_sub_corr.xmg");
      // const djvec_t eff_P5P5_pi0_sub=effective_mass(P5P5_pi0_sub);
      // const djack_t aM_pi0_sub=constant_fit(eff_P5P5_pi0_sub,tmin,TH,"plots/eff_mass_P5P5_pi0_sub.xmg");
      // const djack_t M_rat=aM_pi0_sub/aM;
      // cout<<"MPi0/Mpi+: "<<M_rat.ave_err()<<endl;

	//   {
    // 	jack_fit_t fitter;
    // 	djvec_t pars(3);
    // 	fitter.add_fit_par(pars[0],"C",0.2,0.1);
    // 	fitter.add_fit_par(pars[1],"M",0.14,0.1);
    // 	fitter.add_fit_par(pars[2],"Off",0.04,0.1);
    
    // 	const size_t tmax=TH;
    // 	for(size_t it=5;it<=tmax;it++)
    //   fitter.add_point(P5P5_pi0[it],[t=(double)it,aM](const vector<double>& pars,const size_t iel)
    //   {
    // 	return pars[0]+pars[1]*cosh(pars[2]*(TH-t));
    //   });
    // fitter.fit();
    
    // const djvec_t P5P5_pi0_sub=djvec_t(P5P5_pi0-P5P5_pi0[TH])+P5P5_pi0[TH].ave();
    //   const djvec_t P5P5_pi0_rat=P5P5_pi0/P5P5_pi0.shift(1);
    //   const djvec_t P5P5_pi0_log=log(P5P5_pi0_rat);
    //   P5P5_pi0_log.ave_err().write("plots/P5P5_pi0_corr_log.xmg");
    // cout<<"Pi0: "<<pars.ave_err()<<endl;
    //   }    
 
  if(file_exists("jacks/bubble_S0P5_OS_ingr"))
    {
      const djvec_t S0P5_S0P5_OS=load_bubble("S0P5_OS","S0P5_OS",true,true);
      S0P5_S0P5_OS.ave_err().write("plots/S0P5_S0P5_OS.xmg");
      
      const djvec_t P5P5_OS=read_djvec("jacks/P5P5_OS",T,0).symmetrized();
      P5P5_OS.ave_err().write("plots/P5P5_OS_corr.xmg");
      const djvec_t eff_P5P5_OS=effective_mass(P5P5_OS);
      const djack_t aM_OS=constant_fit(eff_P5P5_OS,tmin,TH,"plots/eff_mass_P5P5_OS.xmg");
      
      const djvec_t P5_B2=load_bubble("P5_B1","P5_B2",true,false);
      P5_B2.ave_err().write("plots/P5_B2.xmg");
      
      const djvec_t P5P5_SS_OS=read_djvec("jacks/P5P5_SS_OS",T,0).symmetrized();
      P5P5_SS_OS.ave_err().write("plots/P5P5_SS_OS.xmg");
      
      // auto fit_OS=bind(fit,aM_OS,eff_P5P5_OS,_1,_2);
      
      // const djvec_t rat_conn_OS=P5P5_SS_OS/P5P5_OS;
      // rat_conn_OS.ave_err().write("plots/ratio_conn_OS.xmg");
      
      // const djvec_t rat_disco_OS=S0P5_S0P5_OS/P5P5_OS;
      // rat_disco_OS.ave_err().write("plots/ratio_disco_OS.xmg");
      
      // const djvec_t rat_combo_OS=(rat_conn_OS-rat_disco_OS)*norm;
      // rat_combo_OS.ave_err().write("plots/ratio_combo_OS.xmg");
      
      // const djvec_t pars_conn_OS=fit_OS(rat_conn_OS,"conn_OS");
      // const djvec_t pars_disco_OS=fit_OS(rat_disco_OS,"disco_OS");
      
      // const djvec_t eff_slope_combo_OS=effective_slope(rat_combo_OS,eff_P5P5_OS,TH);
      
      // const djack_t l7_OS=constant_fit(eff_slope_combo_OS,9,15,"plots/eff_slope_combo_OS.xmg");
      // cout<<"OS l7: "<<l7_OS.ave_err()<<endl;
    }
  
  if(file_exists("jacks/P5P5_0S_TM"))
    {
      const djvec_t P5P5_0S_conn=read_djvec("jacks/P5P5_0S_TM",T,0).symmetrized();
      P5P5_0S_conn.ave_err().write("plots/P5P5_0S_conn_TM.xmg");
      const djvec_t P5P5_0S_disc=load_bubble("S0P5_TM","P5",true,false,false);
      P5P5_0S_disc.ave_err().write("plots/P5P5_0S_disc_TM.xmg");
      //l7 M2
      const djvec_t P5P5_0S_combo=P5P5_0S_conn-P5P5_0S_disc;
      P5P5_0S_combo.ave_err().write("plots/P5P5_0S_combo.xmg");
      
      const djack_t norm_m2=1/zp_fr_zs*sqr(aml)*af_TM/pow(aM_TM,4);
      
      const djvec_t l7_m2_corr=P5P5_0S_combo/P5P5_TM_SL*ZS_TM*norm_m2;
      l7_m2_corr.ave_err().write("plots/l7_m2.xmg");
      
      
      
      // djack_t Z2,M,DZ2_fr_Z2,SL;
  // two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL,P5P5,P5P5_0S,TH,tmin,TH-1,"/dev/null","plots/rat_P5P5_0S.xmg");
  // cout<<"Slope1: "<<SL.ave_err()<<endl;
  // const djack_t A=Z2*exp(-M*TH)/M;
  // const djack_t A1_fr_A=DZ2_fr_Z2+SL*(1.0/M+TH);
  
    }
  
  return 0;
}
