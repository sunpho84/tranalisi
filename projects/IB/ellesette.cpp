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

// #define OS

template <typename T,
	  typename TV>
T fitSlope2Ansatz(const TV& pars,const T&aM,double t)
{
  T th=tanh(aM*(t-TH));
  return pars[0]+pars[1]*sqr(t-TH)+pars[2]*(t-TH)*th;
}

djvec_t load_bubble(const string TAG,const bool useRe,const bool useIm)
{
  cout<<"======== loading bubble "<<TAG<<" =========="<<endl;
  
  ifstream is("jacks/bubble_"+TAG+"_ingr");
  
  vector<complex<double>> d;
  double r,i;
  while(is>>r>>i)
    d.push_back({r,i});
  
  const int nconfs_max=d.size()/T;
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
  
  djvec_t aveBubble(2);
  djack_t aveBubbleSquare;
  // djvec_t avePerTimeRe(T);
  // djvec_t avePerTimeIm(T);
   for(int iconf=0;iconf<nconfs;iconf++)
     {
       const int iclust=iconf/clust_size;
      
       for(int dt=0;dt<(int)T;dt++)
	 {
  // 	  avePerTimeRe[dt][iclust]+=d[dt+T*iconf].real();
  // 	  avePerTimeIm[dt][iclust]+=d[dt+T*iconf].imag();
	   
	   double r=d[dt+T*iconf].real();
	   double i=d[dt+T*iconf].imag();
	   
	   if(useRe)
	     {
	       aveBubbleSquare[iclust]+=r*r;
	       aveBubble[RE][iclust]+=r;
	     }
	   
	   if(useIm)
	     {
	       aveBubbleSquare[iclust]+=i*i;
	       aveBubble[IM][iclust]+=i;
	     }
	 }
     }
   aveBubble.clusterize(clust_size);
   aveBubbleSquare.clusterize(clust_size);
   djack_t aveBubbleSub;
   if(useRe)
     aveBubbleSub+=sqr(aveBubble[RE]);
   if(useIm)
     aveBubbleSub+=sqr(aveBubble[IM]);
   aveBubbleSub-=aveBubbleSquare;
   aveBubbleSub/=(T*(T-1));
   // aveBubble/=T;
   // aveBubbleSquare/=T;
   
   double test=0;
   int n=0;
   for(int iconf=0;iconf<nconfs;iconf++)
     for(int jconf=0;jconf<nconfs;jconf++)
       if(iconf!=jconf)
	 for(int t1=0;t1<(int)T;t1++)
	   for(int t2=0;t2<(int)T;t2++)
	     if(t1!=t2)
	       {
		 test+=
		   (d[t1+T*iconf]*conj(d[t2+T*jconf])).real();
		 n++;
	       }
   // avePerTimeRe.clusterize(clust_size);
  // avePerTimeIm.clusterize(clust_size);
  // aveBubble.clusterize(clust_size);
  // aveBubble/=T;
  // cout<<"Ave of bubble: "<<aveBubble.ave_err()<<endl;
  
  djvec_t bubbleProd(T);
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      const int iclust=iconf/clust_size;
      
      for(int dt=0;dt<(int)T;dt++)
	for(int t1=0;t1<(int)T;t1++)
	  {
	    int t2=(t1+dt)%T;
	    if(useRe)
	      bubbleProd[dt][iclust]+=d[t1+T*iconf].real()*d[t2+T*iconf].real();
	    if(useIm)
	      bubbleProd[dt][iclust]+=d[t1+T*iconf].imag()*d[t2+T*iconf].imag();
	  }
    }
  
  bubbleProd.clusterize(clust_size);
  bubbleProd/=T;
  
   cout.precision(16);
   cout<<"Bubble[1]: "<<bubbleProd[1].ave_err()<<endl;
   cout<<"Bubble[2]: "<<bubbleProd[2].ave_err()<<endl;
   cout<<"Bubble["<<TH-1<<"]: "<<bubbleProd[TH-1].ave_err()<<endl;
   cout<<"Ave: "<<aveBubbleSub.ave_err()<<" biassed: "<<sqr(aveBubble[IM]/T).ave_err()<<" bias: "<<(aveBubbleSquare/T).ave_err()<<endl;
   cout<<"Test: "<<test/n<<endl;
   
   const djack_t y=bubbleProd[2]/aveBubbleSub-1;
   cout<<"Rel diff: "<<y.ave_err()<<endl;
   
   bubbleProd-=
    aveBubbleSub;
   
   bubbleProd-=(djack_t)bubbleProd[TH];
   
   
  // djvec_t aveBubbleProd(T);
  // for(int dt=0;dt<(int)T;dt++)
  //   for(int t1=0;t1<(int)T;t1++)
  //     {
  // 	int t2=(t1+dt)%T;
  // 	//if(RI==0)
  // 	  aveBubbleProd[dt]+=avePerTimeRe[t1]*avePerTimeRe[t2];
  // 	  //else
  // 	  aveBubbleProd[dt]+=avePerTimeIm[t1]*avePerTimeIm[t2];
  //     }
  // aveBubbleProd/=T;
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
  TH=T/2.0;
  aml=input.read<double>("aml");
  zp_fr_zs=input.read<double>("zp_fr_zs");
  
  const djvec_t S0P5_S0P5=load_bubble("S0P5_TM",true,true);
  S0P5_S0P5.ave_err().write("plots/S0P5_S0P5_TM.xmg");
  
  const djvec_t P5P5=read_djvec("jacks/P5P5_TM",T,0).symmetrized();
  P5P5.ave_err().write("plots/P5P5_corr.xmg");
  const djvec_t eff_P5P5=effective_mass(P5P5);
  
  djack_t aM,Z2_P5;
  two_pts_fit(Z2_P5,aM,P5P5,TH,18,TH,"plots/eff_mass_P5P5.xmg");
  const djack_t af=2*aml*sqrt(Z2_P5)/sqr(aM);
  
  cout<<"aM: "<<aM.ave_err()<<endl;
  cout<<"af: "<<af.ave_err()<<endl;
  
  const djack_t norm=-sqr(af)*sqr(aml)/pow(aM,3);
  cout<<"Norm: "<<norm.ave_err()<<endl;
  
  // const djvec_t P5P5_disco=load_bubble_P5P5();
  // P5P5_disco.ave_err().write("plots/P5P5_disco.xmg");
  // const djvec_t eff_P5P5_disco=effective_mass(P5P5_disco);
  // eff_P5P5_disco.ave_err().write("plots/eff_P5P5_disco.xmg");
  
  
  const djvec_t P5P5_SS=read_djvec("jacks/P5P5_SS_TM",T,0).symmetrized();
  P5P5_SS.ave_err().write("plots/P5P5_SS_TM.xmg");
  
  const djvec_t P5P5_0S=read_djvec("jacks/P5P5_0S_TM",T,0).symmetrized();
  P5P5_0S.ave_err().write("plots/P5P5_0S_TM.xmg");
  djack_t Z2,M,DZ2_fr_Z2,SL;
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL,P5P5,P5P5_0S,TH,tmin,TH-1,"/dev/null","plots/rat_P5P5_0S.xmg");
  cout<<"Slope1: "<<SL.ave_err()<<endl;
  const djack_t A=Z2*exp(-M*TH)/M;
  const djack_t A1_fr_A=DZ2_fr_Z2+SL*(1.0/M+TH);
  
  const djvec_t rat_conn=P5P5_SS/P5P5;
  rat_conn.ave_err().write("plots/ratio_conn_TM.xmg");
  
  djvec_t rat_conn_sub=rat_conn;
  for(size_t t=0;t<=TH;t++)
    {
      double dt=(double)TH-t;
      rat_conn_sub[t]-=sqr(SL)*sqr(dt);
    }
  rat_conn_sub.ave_err().write("plots/ratio_conn_sub_TM.xmg");
  const djvec_t eff_slope_conn_sub=effective_slope(rat_conn_sub,eff_P5P5,TH);
  eff_slope_conn_sub.ave_err().write("plots/eff_slope_conn_sub_TM.xmg");
  const djvec_t eff_slope_conn_sub_sub=eff_slope_conn_sub-2*A1_fr_A*SL;
  const djack_t conn_l7_fit=constant_fit(eff_slope_conn_sub_sub,tmin,TH-1,"plots/eff_slope_conn_sub_sub_TM.xmg");
  const djack_t l7=conn_l7_fit*norm;
  cout<<"l7: "<<l7.ave_err()<<endl;
  
  const djvec_t rat_disco=S0P5_S0P5/P5P5;
  rat_disco.ave_err().write("plots/ratio_disco.xmg");
  
  const djvec_t rat_combo=rat_conn-rat_disco/sqr(zp_fr_zs);
  rat_combo.ave_err().write("plots/ratio_combo.xmg");
  
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
  
  const djvec_t pars_conn=fit(aM,eff_P5P5,rat_conn,"conn");
  const djvec_t pars_disco=fit(aM,eff_P5P5,rat_disco,"disco");
  
  const djvec_t pars_combo=fit(aM,eff_P5P5,rat_combo*norm,"combo");
  const djack_t l7combo=pars_combo[2];
  const djack_t l7sub=l7combo-T*pars_combo[1];
  cout<<"l7 combo: "<<l7combo.ave_err()<<endl;
  cout<<"l7 sub: "<<l7sub.ave_err()<<endl;
  
  const djack_t estim_zp_fr_zs=sqrt(pars_disco[1]/pars_conn[1]);
  const djack_t fact=estim_zp_fr_zs/zp_fr_zs;
  cout<<"Estimator: "<<fact.ave_err()<<endl;
  
  const djvec_t rat_eff_combo=norm*(rat_conn-rat_disco/sqr(estim_zp_fr_zs));
  const djvec_t pars_eff_combo=fit(aM,eff_P5P5,rat_eff_combo,"eff_combo");
  const djack_t l7_eff_combo=pars_eff_combo[2];
  cout<<"l7 eff combo: "<<l7_eff_combo.ave_err()<<endl;
  
  const djack_t l7_eff_combo_sub=l7_eff_combo-pars_eff_combo[1]*T;
  cout<<"l7 eff sub: "<<l7_eff_combo_sub.ave_err()<<endl;
  
      const djvec_t P5_B1=load_bubble("P5_B1",false,true);
      P5_B1.ave_err().write("plots/P5_B1.xmg");
      effective_mass(P5_B1).ave_err().write("plots/eff_mass_P5_B1.xmg");
      
  if(file_exists("jacks/bubble_S0P5_OS_ingr"))
    {
      const djvec_t S0P5_S0P5_OS=load_bubble("S0P5_OS",true,true);
      S0P5_S0P5_OS.ave_err().write("plots/S0P5_S0P5_OS.xmg");
      
      const djvec_t P5P5_OS=read_djvec("jacks/P5P5_OS",T,0).symmetrized();
      P5P5_OS.ave_err().write("plots/P5P5_OS_corr.xmg");
      const djvec_t eff_P5P5_OS=effective_mass(P5P5_OS);
      const djack_t aM_OS=constant_fit(eff_P5P5_OS,tmin,TH,"plots/eff_mass_P5P5_OS.xmg");
      
      const djvec_t P5_B2=load_bubble("P5_B1",true,false);
      P5_B2.ave_err().write("plots/P5_B2.xmg");
      
      const djvec_t P5P5_SS_OS=read_djvec("jacks/P5P5_SS_OS",T,0).symmetrized();
      P5P5_SS_OS.ave_err().write("plots/P5P5_SS_OS.xmg");
      
      auto fit_OS=bind(fit,aM_OS,eff_P5P5_OS,_1,_2);
      
      const djvec_t rat_conn_OS=P5P5_SS_OS/P5P5_OS;
      rat_conn_OS.ave_err().write("plots/ratio_conn_OS.xmg");
      
      const djvec_t rat_disco_OS=S0P5_S0P5_OS/P5P5_OS;
      rat_disco_OS.ave_err().write("plots/ratio_disco_OS.xmg");
      
      const djvec_t rat_combo_OS=(rat_conn_OS-rat_disco_OS)*norm;
      rat_combo_OS.ave_err().write("plots/ratio_combo_OS.xmg");
      
      const djvec_t pars_conn_OS=fit_OS(rat_conn_OS,"conn_OS");
      const djvec_t pars_disco_OS=fit_OS(rat_disco_OS,"disco_OS");
      
      const djvec_t eff_slope_combo_OS=effective_slope(rat_combo_OS,eff_P5P5_OS,TH);
      
      const djack_t l7_OS=constant_fit(eff_slope_combo_OS,9,15,"plots/eff_slope_combo_OS.xmg");
      cout<<"OS l7: "<<l7_OS.ave_err()<<endl;
    }
  
  return 0;
}
