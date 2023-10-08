#include <tranalisi.hpp>

#include <loadCorrCov.hpp>

#include <trambacco.hpp>

using Real=double;

const double RE_exp=2;

inline Real thetaSigma(const Real& sigma,
			 const Real& x)
{
  //return 0.5+atan(x/sigma)/M_PI;
  return 1/(1+exp(-x/sigma));
}

inline Real Kl(const Real& sigma,
	       const Real& x)
{
  return sqr(1-x*x)*thetaSigma(sigma,1-x)/x;
}

inline Real Kt(const Real& sigma,
	       const Real& x)
{
  return (1+2*x*x)*Kl(sigma,x);
}

template <size_t I>
vector<pair<djack_t,djack_t>> multiFit(const djvec_t& in,
				       const size_t& T,
				       const string& basePath,
				       const size_t& glbTMin,
				       const size_t& glbTMax)
{
  return {};
}

template <size_t I,
	  typename...Args>
auto multiFit(const djvec_t& in,
	      const size_t& T,
	      const string& basePath,
	      const size_t& glbTMin,
	      const size_t& glbTMax,
	      const size_t& tMin,
	      const size_t& tMax,
	      const string& path,
	      const Args&...args)
{
  const auto [Z2,M]=two_pts_fit(in,T/2,tMin,tMax,basePath+path);
  
  djvec_t corr=in;
  for(size_t iT=0;iT<corr.size();iT++)
    corr[iT]-=two_pts_corr_fun(Z2.ave(),M.ave(),T/2.0,iT,+1);
  
  auto tmp=multiFit<I+1>(corr,T,basePath,glbTMin,glbTMax,args...);
  
  tmp.push_back({Z2,M});
  
  if constexpr(I==0)
    {
      djvec_t p(2*tmp.size());
      
      //parameters to fit
      jack_fit_t fit;
      for(size_t i=0;i<tmp.size();i++)
	{
	  const auto& [Z2Guess,MGuess]=tmp[i];
	  
	  const double Mmin=(MGuess.ave()+((i==tmp.size()-1)?0:tmp[i+1].second.ave()))/2;
	  const double Mmax=(((i==0)?(2*MGuess.ave()):tmp[i-1].second.ave())+MGuess.ave())/2;
	  
	  fit.add_fit_par_limits(p[2*i+0],"Z2G"+to_string(i),Z2Guess.ave_err(),
				 Z2Guess.ave()/10,
				 Z2Guess.ave()*10);
	  fit.add_fit_par_limits(p[2*i+1],"MG"+to_string(i),MGuess.ave_err(),
				 Mmin,
				 Mmax);
				 
	  p[2*i+0]=Z2Guess;
	  p[2*i+1]=MGuess;
	  
	  cout<<"Guess: "<<i<<" Z2="<<Z2Guess.ave_err()<<" , M="<<MGuess.ave_err()<<" in ["<<Mmin<<","<<Mmax<<"]"<<endl;
	}
      
      const auto f=
	[T](const auto& p,
	    const double& t)
	{
	  std::decay_t<decltype(p[0])> s;
	  s=0;
	  for(size_t i=0;i<p.size()/2;i++)
	    s+=two_pts_corr_fun(p[2*i+0],p[2*i+1],T/2.0,t,+1);
	  return s;
	};
      
      for(size_t iT=glbTMin;iT<glbTMax;iT++)
	fit.add_point(//numerical data
		      [o=in[iT]]
		      (const vector<double>& p,
		       const int& iel)
		      {
			return o.ave();
		      },
		      //ansatz
		      [&f,iT]
		      (const vector<double>& fp,
		       const int& iel)
		      {
			return f(fp,iT);
		      },
		      //for covariance/error
		      in[iT].err());
      
      // double initCh2=0;
      // for(size_t iT=glbTMin;iT<glbTMax;iT++)
      // 	{
      // 	  const double a=in[iT].ave();
      // 	  const double t=f(p.ave(),iT);
      // 	  const double e=in[iT].err();
      // 	  const double contr=sqr((a-t)/e);
	  
      // 	  cout<<"contr=(sqr(("<<a<<"-"<<t<<")/"<<e<<")"<<endl;
      // 	  initCh2+=contr;
      // 	}
      // cout<<"Initial ch2: "<<initCh2<<endl;
      
      // fit_debug=true;
      fit.fit();
      
      // double finalCh2=0;
      // for(size_t iT=glbTMin;iT<glbTMax;iT++)
      // 	finalCh2+=sqr((in[iT].ave()-f(p.ave(),iT))/in[iT].err());
      // cout<<"Final ch2: "<<finalCh2<<endl;
      
      grace_file_t plot(basePath+"glbFit.xmg");
      plot.write_vec_ave_err(effective_mass(in).ave_err());
      plot.write_polygon([&](const double& x)
      {
	return effective_mass(f(p,x),f(p,x+1),x,T/2);
      },glbTMin,glbTMax,grace::BLUE);
      
      for(size_t i=0;i<tmp.size();i++)
	{
	  auto& [Z2Out,MOut]=tmp[i];
	  Z2Out=p[2*i+0];
	  MOut=p[2*i+1];
	}
    }
  
  return tmp;
}

int main()
{
  set_njacks(50);
  if constexpr(std::is_same_v<PrecFloat,Real>)
    PrecFloat::setDefaultPrecision(128);
  
  const map<char,Real> zVlist{{'Z',0.697},{'B',0.706379},{'C',0.725404},{'D',0.744108},{'E',0.7663}};
  const map<char,Real> zAlist{{'Z',0.7357},{'B',0.74294},{'C',0.75830},{'D',0.77395},{'E',0.7921}};
  
#if 1
  const size_t T=224;
  const size_t L=112;
  const string ensName="C.06.112";
  const string pref="mix_fixed";
#endif
  
#if 0
  const size_t T=160;
  const size_t L=80;
  const string ensName="C.06.80";
  const string pref="";
#endif
  
#if 0
  const size_t T=128;
  const size_t L=64;
  const string ensName="B.72.64";
  const string pref="mix_fixed";
#endif
  
#if 0
  const size_t T=192;
  const size_t L=96;
  const string ensName="D.54.96";
  const string pref="mix";
#endif
  
  const double z=zAlist.at(ensName[0]);//0.75830;
  
  const string baseOut=
    "/home/francesco/QCD/LAVORI/TAU_VUS/";
  
  auto loadCorrCov=
    [&ensName,&baseOut](const string& chann)
    {
      const string baseIn=
	"/home/francesco/QCD/LAVORI/GM3/";
      
      const auto [readT,corr,cov]=
	::loadCorrCov(baseIn+"/"+ensName+"/data/"+chann,false);
      
      if(T!=readT)
	CRASH("T=%zu, read %zu",T,readT);
      
      corr.ave_err().write(baseOut+ensName+"/"+chann+".xmg");
      effective_mass(corr).ave_err().write(baseOut+ensName+"/eff_mass_"+chann+".xmg");
      
      return make_tuple(corr,cov);
    };
  
  auto loadCorr=
    [&loadCorrCov](const string& chann)
    {
      const auto [corr,cov]=
	loadCorrCov(chann);
      
      return corr;
    };
    
  const djvec_t cP5P5_ll=
    loadCorr(pref + "_l_l_TM_P5P5");
  
  /// Estimates the lattice spacing
  // {
    const djack_t aInv=
      estimateA(aMlist.at('C'),cP5P5_ll,T,baseOut+"/"+ensName+"/pion.xmg");
    cout<<"a^-1: "<<aInv.ave_err()<<" GeV^-1"<<endl;
    cout<<"a: "<<(0.197/aInv).ave_err()<<" fm"<<endl;
  // }
  
  // const djack_t aInFm(gauss_filler_t{0.06820,0.00013,23523523});
  // const djack_t aInv=0.197327/aInFm;
  
#define LOAD(TAG) \
  const auto [c ## TAG ## _ls,cov ## TAG ## _ls]=\
    loadCorrCov(pref + "_l_s1_TM_" #TAG)
  
  LOAD(V0V0);
  LOAD(VKVK);
  LOAD(A0A0);
  LOAD(AKAK);
  LOAD(P5P5);
  
#undef LOAD
  
  const double KaonStarSub=0;
  cout<<"Subtracting the lower state: "<<KaonStarSub<<endl;
  
  const djack_t MPion=constant_fit(effective_mass(cP5P5_ll),30,60,baseOut+ensName+"/pion_fit.xmg");
  const djack_t MKaon=constant_fit(effective_mass(cP5P5_ls),30,60,baseOut+ensName+"/kaon_fit.xmg");
  
  const auto [Z2KaonStar,MKaonStar]=two_pts_fit(cVKVK_ls,T/2,22,35,baseOut+ensName+"/kaonStar_fit.xmg");
  const auto [Z2K1Star,MK1Star]=two_pts_fit(cAKAK_ls,T/2,20,24,baseOut+ensName+"/k1Star_fit.xmg");
  
  const double EAn=MKaonStar.ave()+MKaonStar.err();
  
  djvec_t cVKVK_ls_sub=cVKVK_ls;
  for(size_t iT=0;iT<cVKVK_ls_sub.size();iT++)
    cVKVK_ls_sub[iT]-=two_pts_corr_fun(Z2KaonStar.ave(),EAn,T/2.0,iT,+1)*KaonStarSub;
  
  const double aPmin=2*M_PI/L;
  
  cout<<"MK: "<<smart_print(MKaon*aInv)<<" GeV"<<endl;
  cout<<"MK*: "<<smart_print(MKaonStar*aInv)<<" GeV"<<endl;
  cout<<"MK1*: "<<smart_print(MK1Star*aInv)<<" GeV"<<endl;
  const djack_t EKPi_nonint_FV=sqrt(sqr(MKaon)+sqr(aPmin))+sqrt(sqr(MPion)+sqr(aPmin));
  cout<<"E(K+Pi) non int FVE: "<<smart_print(EKPi_nonint_FV*aInv)<<" GeV"<<endl;
  
   const auto res=
     multiFit<0>(cVKVK_ls,T,baseOut+ensName+"/",10,35,22,35,"kaonStar0_fit.xmg",10,13,"kaonStar_exc_fit.xmg");
  cout<<"K* from multifit: "<<(res[1].second*aInv).ave_err()<<endl;
  
  grace_file_t AKAK_study(baseOut+ensName+"/eff_mass_AKAK_study.xmg");
  AKAK_study.write_vec_ave_err(effective_mass(cAKAK_ls).ave_err());
  AKAK_study.set_legend("AKAK");
  AKAK_study.new_data_set();
  AKAK_study.set_legend("M\\sK\\N+2M\\s\\xp\\0\\N");
  AKAK_study.write_constant_band(0,60,MKaon+2*MPion);
  
  grace_file_t VKVK_study(baseOut+ensName+"/eff_mass_VKVK_study.xmg");
  VKVK_study.write_vec_ave_err(effective_mass(cVKVK_ls).ave_err());
  const djack_t VKVK_ls_mass=constant_fit(effective_mass(cVKVK_ls),20,30);
  VKVK_study.set_legend("VKVK");
  VKVK_study.new_data_set();
  VKVK_study.set_legend("M\\sK\\N+M\\s\\xp\\0\\N");
  VKVK_study.write_constant_band(0,60,MKaon+MPion);
  VKVK_study.set_legend("E(K+Pi) non int FV");
  VKVK_study.write_constant_band(0,60,EKPi_nonint_FV);
  
  const double EMinIntInGeV=0.5;
  const double EMinInt=EMinIntInGeV/aInv.ave();
  const double EMaxInt=6;
  
  const double sigma=0;
  // for(double sigma=0.06;sigma>=0.01;sigma/=1.1)
    {
      //  const double sigma=0.03;
  
      const size_t tMin=1,tMax=43;
  const Basis<Real> basis(T,tMin,tMax,+1);
  const size_t nT=basis.nT;
  
  djvec_t corr(nT);
  
  MatrixX<Real> corrCov(nT,nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      const size_t t=iT+tMin;
      
      corr[iT]=cVKVK_ls[t]-two_pts_corr_fun(Z2KaonStar.ave(),EAn,T/2.0,t,+1)*KaonStarSub;
      for(size_t iS=0;iS<nT;iS++)
	corrCov(iT,iS)=covVKVK_ls[t+(T/2+1)*(iS+tMin)];
    }
  
  corr.ave_err().write(baseOut+"/"+ensName+"/corrToReco.xmg");
  effective_mass(corr,T/2,0).ave_err().write(baseOut+"/"+ensName+"/eff_mass_corrToReco.xmg");
  
  /// Preconditioner of the problem
  vector<Real> preco(nT);
  for(size_t iT=0;iT<basis.nT;iT++)
    preco[iT]=1/corr[iT].ave();
  
  /// Normalization of the correlation function, when assuming that R(E)=E^N
  djvec_t normOfT(nT);
  const double tol=1e-8;
  for(size_t iT=0;iT<nT;iT++)
    {
      const djack_t c=
	corr[iT]-
	integrateUpToInfinity<double>([&basis,
				       iT,
				       p=res[1].first.ave()/(2*res[1].second.ave()),
				       k=res[1].second.ave(),
				       w=0.050/aInv.ave()](const double& E)
      {
	return p*gauss(k,w,E)*basis(iT,E);
      });
				    
      /// Correlator assuming that R(E)=E^N
      const Real cAss=
	gslIntegrateFromTo([iT,&basis](const double& E)
	{
	  return get(pow(E,RE_exp)*basis(iT,E));
	},MKaonStar.ave(),EMaxInt,tol);
      
      normOfT[iT]=c/get(cAss);
    }
  normOfT.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  /// Time at which to normalize the correlator: we choose 2 GeV^-1
  const size_t tNorm=round(2*aInv.ave());
  cout<<"tNorm: "<<tNorm<<endl;
  
  /// Ansatz of the spectral density
  const auto specAns=
    [norm=normOfT[tNorm].ave(),
     p=res[1].first.ave()/(2*res[1].second.ave()),
     k=res[1].second.ave(),
     w=0.050/aInv.ave()](const Real& E)
    {
      return norm*pow(E,RE_exp)+p*gauss(k,w,E);
    };
  
  constexpr double mTau=1.77686;
  
  const djack_t rescGeV=
    sqr(z)*sqr(aInv)*aInv.ave()*24*sqr(M_PI)/cube(mTau);
  
  const auto targetFunction=
    [sigma,amTau=mTau/aInv.ave()](const Real& E)
    {
      const Real x=E/amTau;
      
      return Kt(sigma,x);
    };
  
  const TrambaccoFunctional tf=
    getTrambaccoFunctional<Real>(basis,
				 specAns,
				 targetFunction,
				 EMinInt,EMaxInt,corrCov,preco);
  
  /// Track the results
  djack_t R;
  VectorX<Real> g;
  
  const Real Kcontr=
    Z2KaonStar.ave()/(2*EAn)*targetFunction(EAn)*KaonStarSub;
  
  const auto reconstruct=
    [&](const double& statOverSyst=1)
    {
      jack_t<Real> Rprec;
      
      g=tf.getG(preco,statOverSyst);
      Rprec=0;
      for(size_t iT=0;iT<nT;iT++)
	jackCall([](Real& Rprec,
		    const Real& gIt,
		    const Real& corrIt)
	{
	  Rprec+=gIt*corrIt;
	},Rprec,g[iT],corr[iT]);
      
      Rprec+=Kcontr;
      Rprec*=rescGeV;
      
      jackCall([](double& R,
		  const Real& Rprec)
      {
	R=get(Rprec);
      },R,Rprec);
    };
  
  // Performs the stability check
  grace_file_t stabilityPlot(baseOut+"/"+ensName+"/stab"+to_string(sigma)+".xmg");
  vector<pair<double,djack_t>> stab;
  for(double statOverSyst=512;statOverSyst>=1e-10;statOverSyst/=2)
    {
      reconstruct(statOverSyst);
      stab.emplace_back(log(statOverSyst),R);
    }
  
  // Final reconstruction
  const double lambdaExp=-2.5;
  reconstruct(exp(lambdaExp));
  
  Real KcontrReco=0;
  for(size_t iT=0;iT<nT;iT++)
    KcontrReco+=g[iT]*basis(iT,EAn);
  KcontrReco*=Z2KaonStar.ave()/(2*EAn)*KaonStarSub;
  cout<<" Contr to be added: "<<Kcontr*rescGeV.ave()<<" "<<KcontrReco*rescGeV.ave()<<" "<<Kcontr/KcontrReco-1<<endl;
  
  for(const auto& [x,y] : stab)
    stabilityPlot.write_ave_err(x,y.ave_err());
  stabilityPlot.write_constant_band(lambdaExp-0.1,lambdaExp+0.1,R);
  
  /// Print coefficients
  grace_file_t gPlot(baseOut+"/"+ensName+"/g.xmg");
  for(size_t iT=0;iT<nT;iT++)
    gPlot.write_xy(iT,get(g[iT]));
  
  const auto interpReco=
    [&g,&aInv,&nT,&basis](const Real& EInGeV)
    {
      const Real E=EInGeV/aInv.ave();
      Real s=0;
      for(size_t iT=0;iT<nT;iT++)
	s+=g[iT]*basis(iT,E);
      
      return s;
    };
  
  // Plots the reconstruction of the target function
  grace_file_t recoPlot(baseOut+"/"+ensName+"/reco.xmg");
  recoPlot.write_line([&](const Real& x)
  {
    return get(interpReco(x)*x);
  },EMinInt*aInv.ave(),EMaxInt*aInv.ave(),1000);
  
  recoPlot.write_line([&](const Real& x)
  {
    const Real EInGeV=x;
    const Real E=EInGeV/aInv.ave();
    
    return get(targetFunction(E)*x);
  },EMinInt*aInv.ave(),EMaxInt*aInv.ave(),1000);
  
  // Plots the relative erorr of the reconstruction of the target function
  grace_file_t errPlot(baseOut+"/"+ensName+"/err.xmg");
  errPlot.write_line([&](const Real& x)
  {
    const double reco=get(interpReco(x)*x);
    
    const Real EInGeV=x;
    const Real E=EInGeV/aInv.ave();
    
    const double targ=get(targetFunction(E)*x);
    
    return reco-targ;
  },EMinInt*aInv.ave(),EMaxInt*aInv.ave(),1000);
  
  cout<<sigma<<" "<<R.ave_err()<<endl;
    }
  
  return 0;
}
