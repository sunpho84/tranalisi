#include <tranalisi.hpp>

const size_t T=96;
const size_t ntot_cols=2;
range_t confs_range;
double am;

djvec_t load(const string& tag,const size_t offset,const size_t col,const size_t par)
{
  const djvec_t data=read_conf_set_t("out/%04d/mes_contr_"+tag,confs_range,ntot_cols,{0,1},T,VERBOSE);
  
  const size_t i=col+2*offset;
  
  return data.subset(i*T,(i+1)*T-1).symmetrized(par);
}

int main()
{
  set_njacks(10);
  
  input_file_t input("input.txt");
  am=input.read<double>("am");
  const double kappa=input.read<double>("Kappa");
  input.expect("ConfRange");
  confs_range.start=input.read<size_t>();
  confs_range.each=input.read<size_t>();
  confs_range.end=input.read<size_t>();
  const size_t RE=0,IM=1;
  const size_t P5P5=0,V0P5=1;
  const size_t EVN=1,ODD=-1;
  const djvec_t cP5P5=load("NOINS",P5P5,RE,EVN);
  const djvec_t cV0P5=load("NOINS",V0P5,IM,ODD);
  const djvec_t cP5PP5=2*load("PINS",P5P5,IM,EVN); //- from tau3, -from i
  const djvec_t cV0PP5=-2*load("PINS",V0P5,RE,ODD);
  const djvec_t cP5SP5=2*load("SINS",P5P5,RE,EVN);
  const djvec_t cV0SP5=2*load("SINS",V0P5,IM,ODD);
  
  cP5P5.ave_err().write("plots/cP5P5.xmg");
  const size_t tmin=19,tmax=20;
  // const djvec_t eff_mass_P5=effective_mass(cP5P5);
  // eff_mass_P5.ave_err().write("plots/meff_P5.xmg");
  // const djvec_t eff_slope_P=effective_slope(cP5PP5/cP5P5,eff_mass_P5,T/2);
  // const djvec_t eff_slope_S=effective_slope(cP5SP5/cP5P5,eff_mass_P5,T/2);
  djack_t Z2,aMPi,DPZ2_fr_Z2,SL_P,DSZ2_fr_Z2,SL_S;
  two_pts_with_ins_ratio_fit(Z2,aMPi,DPZ2_fr_Z2,SL_P,cP5P5,cP5PP5,T/2,tmin,tmax,"plots/Meff_P5_Pfit.xmg","plots/seff_P5_Pfit.xmg",+1);
  two_pts_with_ins_ratio_fit(Z2,aMPi,DSZ2_fr_Z2,SL_S,cP5P5,cP5SP5,T/2,tmin,tmax,"plots/Meff_P5_Sfit.xmg","plots/seff_P5_Sfit.xmg",+1);
  
  cout<<"SL_S: "<<SL_S.ave_err()<<endl;
  cout<<"SL_P: "<<SL_P.ave_err()<<endl;
  const djack_t aF_Pi=2*am*sqrt(Z2)/sqr(aMPi);
  const djack_t dSaM_Pi_fr_aMPi=SL_S/aMPi;
  const djack_t dPaM_Pi_fr_aMPi=SL_P/aMPi;
  const djack_t dSaF_Pi_fr_aFPi=1/am+DPZ2_fr_Z2/2-2*dSaM_Pi_fr_aMPi;
  const djack_t dPaF_Pi_fr_aFPi=DPZ2_fr_Z2/2-2*dPaM_Pi_fr_aMPi;
  cout<<"dSaF_Pi/aF_Pi: "<<dSaF_Pi_fr_aFPi.ave_err()<<endl;
  cout<<"dPaF_Pi/aF_Pi: "<<dPaF_Pi_fr_aFPi.ave_err()<<endl;
  const djack_t dS_X_fr_X=dSaM_Pi_fr_aMPi-dSaF_Pi_fr_aFPi;
  const djack_t dP_X_fr_X=dPaM_Pi_fr_aMPi-dPaF_Pi_fr_aFPi;
  
  // cV0P5.ave_err().write("plots/cV0P5.xmg");
  
  // cP5PP5.ave_err().write("plots/cP5PP5.xmg");
  // cV0PP5.ave_err().write("plots/cV0PP5.xmg");
  
  // const djvec_t ratV0P5=cV0PP5/cV0P5;
  
  // Correction
  
  // djvec_t ratV0P5_adj=ratV0P5;
  // for(size_t t=0;t<T/2;t++)
  //   {
  //     const djack_t m=eff_mass_P5[t];
  //     const djack_t a=exp(-m*t),b=exp(-m*(T-t));
  //     const djack_t corr=(t*a-(T-t)*b)/(a-b);
      
  //     ratV0P5_adj[t]+=eff_slope_P[t]*corr;
  //   }
  
  // grace_file_t plot("plots/ratV0P5.xmg");
  // plot.set_no_line();
  // const djack_t dm=1/constant_fit(ratV0P5,tmin,tmax);
  
  //plot.write_constant_band(tmin,tmax,1/dm,grace::RED);
  // plot.write_vec_ave_err(ratV0P5.ave_err(),grace::RED,grace::SQUARE);
  // plot.write_constant_band(tmin,tmax,1/dm_adj,grace::BLACK);
  // plot.write_constant_band(0,tmin,1/dm_adj,grace::MAROON);
  // plot.set_no_line();
  // plot.write_vec_ave_err(ratV0P5_adj.ave_err(),grace::VIOLET,grace::SQUARE);
  
  // const djvec_t dm_lin_pars=poly_fit(ratV0P5,1,tmin,tmax);
  // const djack_t dm_lin=1/dm_lin_pars[0];
  // const auto poly=[&dm_lin_pars](const double x)
  // {
  //   return poly_eval(dm_lin_pars,x);
  // };
  // plot.write_polygon(poly,tmin,tmax,grace::GREEN4);
  // plot.write_polygon(poly,0,tmin,grace::GREEN);
  
  // const djack_t kappa_prime=kappa-2*dm*sqr(kappa);
  // const djack_t kappa_prime_lin=kappa-2*dm_lin*sqr(kappa);
  // const djack_t kappa_prime_adj=kappa-2*dm_adj*sqr(kappa);

  // cout<<"dm adj: "<<smart_print(dm_adj.ave_err())<<endl;
  cout<<"kappa: "<<smart_print(kappa)<<endl;
  // cout<<"kappa prime: "<<smart_print(kappa_prime.ave_err())<<endl;
  // cout<<"kappa prime lin: "<<smart_print(kappa_prime_lin.ave_err())<<endl;
  // cout<<"kappa prime adj: "<</*smart_print*/(kappa_prime_adj.ave_err())<<endl;
  // const djack_t dm_induced=1/(2*kappa_prime_adj)-1/(2*kappa);
  // cout<<"dm_induced: "<<smart_print(dm_induced.ave_err())<<endl;
  
  const djvec_t mPCAC_correl=symmetric_derivative(cV0P5)/(2*cP5P5);
  const djack_t mPCAC=constant_fit(mPCAC_correl,tmin,tmax,"plots/mPCAC.xmg");
  const djvec_t dPmPCAC_fr_mPCAC_correl=symmetric_derivative(cV0PP5)/symmetric_derivative(cV0P5)-cP5PP5/cP5P5;
  const djvec_t dSmPCAC_fr_mPCAC_correl=symmetric_derivative(cV0SP5)/symmetric_derivative(cV0P5)-cP5SP5/cP5P5;
  const djack_t dPmPCAC_fr_mPCAC=constant_fit(dPmPCAC_fr_mPCAC_correl,tmin,tmax,"plots/dPmPCAC_fr_mPCAC.xmg");
  const djack_t mPCAC_fr_dPmPCAC=constant_fit(1/dPmPCAC_fr_mPCAC_correl,tmin,tmax,"plots/mPCAC_fr_dPmPCAC.xmg");
  const djack_t dSmPCAC_fr_mPCAC=constant_fit(dSmPCAC_fr_mPCAC_correl,tmin,tmax,"plots/dSmPCAC_fr_mPCAC.xmg");
  cout<<"mPCAC: "<<mPCAC.ave_err()<<endl;
  cout<<"dPmPCAC/mPCAC: "<<dPmPCAC_fr_mPCAC.ave_err()<<endl;
  cout<<"1/(mPCAC/dPmPCAC): "<<(1/mPCAC_fr_dPmPCAC).ave_err()<<endl;
  cout<<"mPCAC/dPmPCAC: "<<mPCAC_fr_dPmPCAC.ave_err()<<endl;
  cout<<"dSmPCAC/mPCAC: "<<dSmPCAC_fr_mPCAC.ave_err()<<endl;
  
  const djack_t dPm_adj=-mPCAC_fr_dPmPCAC;
  const djack_t kappa_prime_adj=kappa-2*dPm_adj*sqr(kappa);
  cout<<"kappa prime adj: "<</*smart_print*/(kappa_prime_adj.ave_err())<<endl;
  
  
  // const djack_t dPm_adj=-mPCAC/dPmPCAC;
  // cout<<"dPm from mPCAC: "<<dPm_adj.ave_err()<<endl;
  // const djack_t kappa_prime_adj2=kappa-2*dPm_adj*sqr(kappa);
  
  // grace_file_t meff_P5_adj("plots/meff_P5_adj.xmg");
  // const djvec_t cP5P5_adj2=cP5P5+dPm_adj*cP5PP5;
  // meff_P5_adj.write_vec_ave_err(effective_mass(cP5P5).ave_err(),grace::BLACK,grace::CIRCLE);
  // meff_P5_adj.write_vec_ave_err(effective_mass(cP5P5_adj2).ave_err(),grace::RED,grace::SQUARE);
  
  // const double rP5_exp=(135.0/130.4);
  // const double a=0.091/0.197;
  // djack_t aM_P5,z2_P5;
  // two_pts_fit(z2_P5,aM_P5,cP5P5,T/2,tmin,tmax,"plots/aM_P5.xmg");
  // const djack_t M_P5=aM_P5/a;
  // const djack_t aF_Pi=2*am*sqrt(z2_P5)/(aM_P5*sinh(aM_P5));
  // const djack_t F_Pi=aF_Pi/a;
  // const djack_t rP5=M_P5/F_Pi;

  // djack_t aM_P5_adj,z2_P5_adj;
  // two_pts_fit(z2_P5_adj,aM_P5_adj,cP5P5_adj2,T/2,tmin,tmax,"plots/aM_P5_adj.xmg");
  // const djack_t M_P5_adj=aM_P5_adj/a;
  // const djack_t F_Pi_adj=aF_Pi_adj/a;
  // const djack_t rP5_adj=M_P5_adj/F_Pi_adj;
  // cout<<" M_P5: "<<smart_print(M_P5.ave_err())<<" GeV"<<endl;
  // cout<<" M_P5_adj: "<<smart_print(M_P5_adj.ave_err())<<" GeV"<<endl;
  // cout<<" fPi: "<<F_Pi.ave_err()<<" GeV"<<endl;
  // cout<<" fPi adj: "<<F_Pi_adj.ave_err()<<" GeV"<<endl;
  // cout<<" rP5: "<<rP5.ave_err()<<endl;
  // cout<<" rP5_adj: "<<smart_print(rP5_adj.ave_err())<<endl;
  // cout<<" rP5_exp: "<<rP5_exp<<endl;
  
  // const djack_t am_corr=am*sqr(rP5_exp/rP5_adj);
  // cout<<" am_phys_line: "<<smart_print(am_corr.ave_err())<<endl;
  
  if(file_exists("conndisco.txt"))
    {
      const int nConfs=28;
      
      ifstream conndisco("conndisco.txt");
      vector<double> conn(nConfs*T),disco(nConfs);
      const size_t V=48*48*48*96;
      for(size_t iConf=0;iConf<nConfs;iConf++)
	{
	  if(not (conndisco>>disco[iConf])) CRASH("reading disco");
	  disco[iConf]*=V;
	  for(size_t t=0;t<T;t++)
	    if(not (conndisco>>conn[t+T*iConf])) CRASH("reading conn");
	}
      
      djack_t d;
      djvec_t c(T);
      djvec_t dc(T);
      jackknivesFill(nConfs,[&dc,&disco,&conn,&c,&d](const size_t& iConf,const size_t& iClust,const double& w)
      {
	d[iClust]+=disco[iConf]*w;
	for(size_t t=0;t<T;t++)
	  {
	    dc[t][iClust]+=disco[iConf]*conn[t+T*iConf]*w;
	    c[t][iClust]+=conn[t+T*iConf]*w;
	  }
      });
      
      const double discoClustSize=(double)nConfs/njacks;
      c.clusterize(discoClustSize);
      dc.clusterize(discoClustSize);
      d.clusterize(discoClustSize);
      
      dc-=d*c;
      c.symmetrize(0);
      dc.symmetrize(0);
      c.ave_err().write("plots/conn.xmg");
      dc.ave_err().write("plots/disco.xmg");
      
      const djvec_t dmPCAC_correl_conn=-symmetric_derivative(cV0PP5)/(2*cP5P5);
      const djvec_t dmPCAC_correl_disco=-symmetric_derivative(dc)/(2*cP5P5);
      dmPCAC_correl_conn.ave_err().write("plots/dmPCAC_correl_conn.xmg");
      dmPCAC_correl_disco.ave_err().write("plots/dmPCAC_correl_disco.xmg");
    }
  
  return 0;
}
