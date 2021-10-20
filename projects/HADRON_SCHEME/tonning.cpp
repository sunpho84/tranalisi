#include <set>

#include <tranalisi.hpp>
#include <unistd.h>

template <typename T>
struct set_ordered_by_insertion
{
  set<T> uniqueList;
  
  vector<T> orderedList;
  
  bool insert(const T& i)
  {
    if(uniqueList.insert(i).second)
      {
	orderedList.push_back(i);
	
	return true;
      }
    else
      return false;
  }
  
  const auto begin() const
  {
    return
      orderedList.begin();
  }
  
  auto end() const
  {
    return
      orderedList.end();
  }
};

tuple<djack_t,djack_t> retune(const djack_t& dM_dmu,
			      const djack_t& dM_dm,
			      const djack_t& dMPCAC_dmu,
			      const djack_t& dMPCAC_dm,
			      const djack_t& dM,
			      const djack_t& dMPCAC)
{
  djack_t dmu,dm;
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      const vector<double> A
	{dM_dmu  [ijack],  dM_dm[ijack],
	 dMPCAC_dmu[ijack],dMPCAC_dm[ijack]};
      
      const vector<double> b{dM[ijack],dMPCAC[ijack]};
      
      const vector<double> x=
	lin_solve<vector<double>,double>(A,b);
      
      dmu[ijack]=x[0];
      dm[ijack]=x[1];
    }
  
  return
    {dmu,dm};
}

int main(int argc, char** argv)
{
  if(argc!=2)
    CRASH("usage: %s <T>",argv[0]);
  
  const size_t T=atoi(argv[1]);
  
  raw_file_t in("tm_tuning","r");
  char line[1024];
  
  using Key=std::array<double,3>;
  
  set_ordered_by_insertion<Key> parsList;
  
  std::map<Key,vector<std::array<double,8>>> raw;
  bool goon=true;
  do
    {
      in.get_line(line);
      istringstream is(line);
      string dum;
      
      size_t conf,flav;
      double mass,kappa,cSW;
      is>>
	dum>>dum>>conf>>
	dum>>dum>>dum>>flav>>
	dum>>dum>>dum>>mass>>
	dum>>dum>>dum>>kappa>>
	dum>>dum>>dum>>cSW;
      
      goon=(not in.feof());
      
      dcomplex l[8];
      
      if(goon)
	{
	  for(size_t t=0;t<T;t++)
	    {
	      in.expect(to_string(t).c_str());
	      
	      for(int i=0;i<8;i++)
		{
		  const double re=in.read<double>();
		  in.expect(",");
		  const double im=in.read<double>();
		  
		  l[i]={re,im};
		}
	      
	      const double P5P5=l[0].real();
	      const double P5P5_OS=l[1].real();
	      const double S0P5_OS=l[2].imag();
	      const double V0P5=l[3].imag();
	      const double P5P5_ins_S=-2*l[4].real();
	      const double V0P5_ins_S=-2*l[5].imag();
	      const double P5P5_ins_P=2*l[6].imag();
	      const double V0P5_ins_P=-2*l[7].real();
	      
	      const Key key={kappa,mass,cSW};
	      
	      parsList.insert(key);
	      
	      raw[key].push_back({P5P5,P5P5_OS,S0P5_OS,V0P5,P5P5_ins_S,V0P5_ins_S,P5P5_ins_P,V0P5_ins_P});
	    }
	  
	  in.get_line(line);
	  in.get_line(line);
	}
    }
  while(goon);
  
  cout<<raw.size()<<endl;
  
  std::vector<double> corr_vec;
  goon=true;
  
  grace_file_t P5P5_plot("P5P5.xmg");
  grace_file_t P5P5_OS_plot("P5P5_OS.xmg");
  grace_file_t mPCACCorr_plot("mPCACCorr.xmg");
  grace_file_t dMPCAC_dmw_plot("dMPCAC_dmw.xmg");
  grace_file_t dMPCAC_dmu_plot("dMPCAC_dmu.xmg");
  grace_file_t dMPi_dmw_plot("dMPi_dmw.xmg");
  grace_file_t dMPi_dmu_plot("dMPi_dmu.xmg");
  grace_file_t V0P5_plot("V0P5.xmg");
  grace_file_t P5P5_der_S_plot("P5P5_der_S.xmg");
  grace_file_t V0P5_der_S_plot("V0P5_der_S.xmg");
  grace_file_t P5P5_der_P_plot("P5P5_der_P.xmg");
  grace_file_t V0P5_der_P_plot("V0P5_der_P.xmg");
  
  grace_file_t P5P5_ins_P_plot("P5P5_ins_P.xmg");
  grace_file_t P5P5_ins_S_plot("P5P5_ins_S.xmg");
  grace_file_t V0P5_ins_P_plot("V0P5_ins_P.xmg");
  grace_file_t V0P5_ins_S_plot("V0P5_ins_S.xmg");
  grace_file_t mPCAC_fun_kappa_plot("mPCAC_fun_kappa.xmg");
  grace_file_t mpi_fun_kappa_plot("mpi_fun_kappa.xmg");
  grace_file_t pion_mass_diff_file("pionMassDiff.xmg");
  grace_file_t pion_mass_file("pionMass.xmg");
  grace_file_t pion_OS_mass_file("pionOSMass.xmg");
  // grace_file_t kprime_file("kprime_corr.xmg");
  
  for(auto& p : parsList)
    {
      const auto& s=raw[p];
      
      const double kappa=p[0];
      const double mw=0.5/kappa;
      const double mu=p[1];
      const double cSW=p[2];
      const size_t n=s.size()/T;
      const size_t dump=0;
      njacks=n-dump;
      
      cout<<"n: "<<njacks<<" , k: "<<kappa<<" , mu: "<<mu<<" , cSW: "<<cSW<<endl;
      cout<<"================================"<<endl;
      
      djvec_t P5P5(T);
      djvec_t P5P5_OS(T);
      djvec_t S0P5_OS(T);
      djvec_t V0P5(T);
      djvec_t P5P5_ins_P(T);
      djvec_t V0P5_ins_P(T);
      djvec_t P5P5_ins_S(T);
      djvec_t V0P5_ins_S(T);
      
      for(size_t iconf=dump;iconf<n;iconf++)
	for(size_t t=0;t<T;t++)
	  {
	    const auto& a=s[t+T*iconf];
	    
	    P5P5[t][iconf-dump]=a[0];
	    P5P5_OS[t][iconf-dump]=a[1];
	    S0P5_OS[t][iconf-dump]=a[2];
	    V0P5[t][iconf-dump]=a[3];
	    P5P5_ins_S[t][iconf-dump]=a[4];
	    V0P5_ins_S[t][iconf-dump]=a[5];
	    P5P5_ins_P[t][iconf-dump]=a[6];
	    V0P5_ins_P[t][iconf-dump]=a[7];
	  }
      
      P5P5.clusterize().symmetrize();
      V0P5.clusterize().symmetrize(-1);
      P5P5_ins_S.clusterize().symmetrize();
      V0P5_ins_S.clusterize().symmetrize(-1);
      P5P5_ins_P.clusterize().symmetrize();
      V0P5_ins_P.clusterize().symmetrize(-1);
      P5P5_OS.clusterize().symmetrize();
      S0P5_OS.clusterize().symmetrize();
      
      P5P5_ins_S_plot.write_vec_ave_err(P5P5_ins_S.ave_err());
      P5P5_ins_P_plot.write_vec_ave_err(P5P5_ins_P.ave_err());
      V0P5_ins_S_plot.write_vec_ave_err(V0P5_ins_S.ave_err());
      V0P5_ins_P_plot.write_vec_ave_err(V0P5_ins_P.ave_err());
      
      const size_t tmin=12,tmax=22;
      
      /////////////////////////////////////////////////////////////////
      
      const djvec_t mPCACCorr=symmetric_derivative(V0P5)/(2*P5P5);
      const djack_t mPCAC=constant_fit(mPCACCorr,tmin,tmax);
      mPCACCorr_plot.write_vec_ave_err(mPCACCorr.ave_err());
      
      cout<<"mPCAC: "<<mPCAC.ave_err()<<endl;
      
      /////////////////////////////////////////////////////////////////
      
      djack_t Z2,aMPi;
      two_pts_fit(Z2,aMPi,P5P5,T/2,tmin,tmax);
      pion_mass_file.write_vec_ave_err(effective_mass(P5P5).ave_err());
      pion_mass_file.write_constant_band(tmin,tmax,aMPi);
      
      const djvec_t eff_pi=effective_mass(P5P5);
      P5P5_plot.write_vec_ave_err(P5P5.ave_err());
      const djvec_t eff_pi_OS=effective_mass(P5P5_OS);
      P5P5_OS_plot.write_vec_ave_err(P5P5_OS.ave_err());
      
      const djvec_t eff_pi_diff=eff_pi_OS-eff_pi;
      pion_mass_diff_file.write_vec_ave_err(eff_pi_diff.ave_err());
      
      mPCAC_fun_kappa_plot.write_ave_err(kappa,mPCAC.ave_err());
      mpi_fun_kappa_plot.write_ave_err(kappa,aMPi.ave_err());
      
      V0P5_plot.write_vec_ave_err(V0P5.ave_err());
      
      const djvec_t P5P5_der_P=P5P5_ins_P/P5P5;
      P5P5_der_P_plot.write_vec_ave_err(P5P5_der_P.ave_err());
      
      const djvec_t P5P5_der_S=P5P5_ins_S/P5P5;
      P5P5_der_S_plot.write_vec_ave_err(P5P5_der_S.ave_err());
      
      const djvec_t V0P5_der_P=V0P5_ins_P/V0P5;
      V0P5_der_P_plot.write_vec_ave_err(V0P5_der_P.ave_err());
      
      const djvec_t V0P5_der_S=V0P5_ins_S/V0P5;
      V0P5_der_S_plot.write_vec_ave_err(V0P5_der_S.ave_err());
      
      /////////////////////////////////////////////////////////////////
      
      const djvec_t dMPCAC_dmw_corr=forward_derivative(V0P5_ins_P)/(2*P5P5)-forward_derivative(V0P5)*P5P5_ins_P/(2*sqr(P5P5));
      const djack_t dMPCAC_dmw=constant_fit(dMPCAC_dmw_corr,tmin,tmax);
      dMPCAC_dmw_plot.write_vec_ave_err(dMPCAC_dmw_corr.ave_err());
      dMPCAC_dmw_plot.write_constant_band(tmin,tmax,dMPCAC_dmw);
      
      const djvec_t dMPCAC_dmu_corr=forward_derivative(V0P5_ins_S)/(2*P5P5)-forward_derivative(V0P5)*P5P5_ins_S/(2*sqr(P5P5));
      const djack_t dMPCAC_dmu=constant_fit(dMPCAC_dmu_corr,tmin,tmax);
      dMPCAC_dmu_plot.write_vec_ave_err(dMPCAC_dmu_corr.ave_err());
      dMPCAC_dmu_plot.write_constant_band(tmin,tmax,dMPCAC_dmu);
      
      const djvec_t dMPi_dmw_corr=-effective_slope(P5P5_ins_P/P5P5,effective_mass(P5P5),T/2);
      const djack_t dMPi_dmw=constant_fit(dMPi_dmw_corr,tmin,tmax);
      dMPi_dmw_plot.write_vec_ave_err(dMPi_dmw_corr.ave_err());
      dMPi_dmw_plot.write_constant_band(tmin,tmax,dMPi_dmw);
      
      const djvec_t dMPi_dmu_corr=-effective_slope(P5P5_ins_S/P5P5,effective_mass(P5P5),T/2);
      const djack_t dMPi_dmu=constant_fit(dMPi_dmu_corr,tmin,tmax);
      dMPi_dmu_plot.write_vec_ave_err(dMPi_dmu_corr.ave_err());
      dMPi_dmu_plot.write_constant_band(tmin,tmax,dMPi_dmu);
      
      /////////////////////////////////////////////////////////////////
      
      const double aMPiTarget=0.105354;
      const djack_t dM=aMPiTarget-aMPi;
      const djack_t d2M=sqr(aMPiTarget)-sqr(aMPi);
      const djack_t dMPCAC=-mPCAC;
      cout<<"mPiTarget: "<<aMPiTarget<<endl;
      // const auto [dmu,dmw]=retune(2*dMPi_dmu*aMPi,2*dMPi_dmw*aMPi,
      // 				  dMPCAC_dmu,dMPCAC_dmw,
      // 				  d2M,dMPCAC);
      const auto [dmu,dmw]=retune(dMPi_dmu,dMPi_dmw,
				  dMPCAC_dmu,dMPCAC_dmw,
				  dM,dMPCAC);
      const djack_t muRetuned=mu+dmu;
      const djack_t kappaRetuned=0.5/(mw+dmw);
      // cout<<" dmu: "<<dmu.ave_err()<<endl;
      // cout<<" dmw: "<<dmw.ave_err()<<endl;
      cout<<"Retuned mu: "<<muRetuned.ave_err()<<endl;
      cout<<"Retuned kappa:  "<<kappaRetuned.ave_err()<<endl;
      cout<<"Retuned kappa ignoring m:  "<<(0.5/(mw+dMPCAC/dMPCAC_dmw)).ave_err()<<endl;
      // const djack_t e=dMPCAC/dMPCAC_dmw;
      // const djack_t f=0.5/(mw+e);
      // const djack_t g=dM/dMPi_dmw;
      // const djack_t h=0.5/(mw+g);
      // cout<<"Retuned kappa fomr mpc:  "<<f.ave_err()<<endl;
      // cout<<"Retuned kappa from m2pi: "<<h.ave_err()<<endl;
      /////////////////////////////////////////////////////////////////
      
      // const djack_t dmwAlt=-mPCAC/dMPCAC_dmw;
      // const djvec_t dmw_corr=-mcr/dmcr_dmw_corr;
      // const djack_t kappa_prime=0.5/(mw+dmwAlt);
      // const djvec_t kappa_prime_corr=0.5/(mw+dmw_corr);
      // cout<<"mcr: "<<mcr.ave_err()<<endl;
      // cout<<"dmAlt: "<<dmwAlt.ave_err()<<endl;
      // cout<<"kappa_prime: "<<kappa_prime.ave_err()<<endl;
      // kprime_file.write_vec_ave_err(kappa_prime_corr.ave_err());
      
      djack_t Z2_OS,aMPi_OS;
      two_pts_fit(Z2_OS,aMPi_OS,P5P5_OS,T/2,tmin,tmax);
      pion_OS_mass_file.write_vec_ave_err(effective_mass(P5P5_OS).ave_err());
      pion_OS_mass_file.write_constant_band(tmin,tmax,aMPi_OS);
      
      const djack_t OS_fr_TM=aMPi_OS/aMPi-1;
      cout<<"Gap: "<<OS_fr_TM.ave_err()<<endl;
      
      const djack_t aFPi=
	2*mu*sqrt(Z2)/(aMPi*sinh(aMPi));
      
      djack_t ZP5,ZS0,aMPiBis;
      two_pts_SL_fit(ZP5,ZS0,aMPiBis,S0P5_OS,P5P5_OS,T/2,tmin,tmax,"S0.xmg",+1,+1);
      if(isnan(ZS0.ave()))
	{
	  djvec_t u=S0P5_OS;
	  for(size_t t=0;t<u.size();t++)
	    u[t]/=two_pts_corr_fun(sqrt(Z2_OS),aMPi_OS,T/2.0,t,+1);
	  ZS0=constant_fit(u,tmin,tmax,"singleFitZS0.xmg");
	}
      const djack_t ZPS0=sqrt(Z2_OS)*ZS0;
      cout<<"ZPS0: "<<ZPS0.ave_err()<<endl;
      
      // const double a=
      // 	0.414;
      
      // const djack_t mPi=
      // 	aMPi/a;
      
      // const djack_t fPi=
      // 	aFPi/a;
      
      const djack_t r=
	aMPi/aFPi;
      
      const double fPiPhys=
	0.1304;
      
      const double mPiPhys=
	0.13498;
      
      const double rPhys=
	mPiPhys/fPiPhys;
      
      const djack_t rDev=
	r/rPhys;
      
      cout<<"aMPi: "<<aMPi.ave_err()<<endl;
      cout<<"aMPi_OS: "<<aMPi_OS.ave_err()<<endl;
      cout<<"aFPi: "<<aFPi.ave_err()<<endl;
      cout<<"r/rPhys: "<<rDev.ave_err()<<endl;
      
      /////////////////////////////////////////////////////////////////
      
      const djvec_t V0P5Prime=
	V0P5+dmw*V0P5_ins_P+dmu*V0P5_ins_S;
      
      const djvec_t P5P5Prime=
	P5P5+dmw*P5P5_ins_P+dmu*P5P5_ins_S;
      
      const djvec_t mPCACPrimeEff=
	symmetric_derivative(V0P5Prime)/(2*P5P5Prime);
      
      cout<<"mPCACPrime: "<<constant_fit(mPCACPrimeEff,tmin,tmax).ave_err()<<endl;
      {
	djack_t Z2Prime,aMPiPrime;
	two_pts_fit(Z2Prime,aMPiPrime,P5P5Prime,T/2,tmin,tmax);
	cout<<"mPiPrime: "<<aMPiPrime.ave_err()<<endl;
	cout<<"mPiPrime: "<<(aMPi+dmu*dMPi_dmu+dmw*dMPi_dmw).ave_err()<<endl;
      }
      
      for(auto& plot : {&P5P5_plot,&mPCACCorr_plot,&dMPCAC_dmw_plot,&dMPCAC_dmu_plot,&dMPi_dmw_plot,&dMPi_dmu_plot,&V0P5_plot,
			&P5P5_der_P_plot,&V0P5_der_P_plot,&P5P5_der_S_plot,&V0P5_der_S_plot,
			&P5P5_ins_P_plot,&V0P5_ins_P_plot,&P5P5_ins_S_plot,&V0P5_ins_S_plot,
			&pion_mass_file,&pion_mass_diff_file,&pion_mass_file,&pion_OS_mass_file})
	{
	  plot->set_settype(grace::XYDY);
	  plot->set_legend(to_string(kappa));
	}
      
      cout<<endl;
    }
  
  return 0;
}
