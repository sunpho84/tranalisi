#include <filesystem>
#include <tranalisi.hpp>

#include <loadCorrCov.hpp>

#include <trambacco.hpp>

using namespace filesystem;

using Real=double;

const double RE_exp=2;

const string baseOut=
  "/home/francesco/QCD/LAVORI/TAU_VUS/";

constexpr double mTau=1.77686;
  
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
	      const path& basePath,
	      const size_t& glbTMin,
	      const size_t& glbTMax,
	      const size_t& tMin,
	      const size_t& tMax,
	      const string& path,
	      const Args&...args)
{
  const auto [Z2,M]=
    two_pts_fit(in,T/2,tMin,tMax,basePath/path);
  
  djvec_t corr=in;
  for(size_t iT=0;iT<corr.size();iT++)
    corr[iT]-=two_pts_corr_fun(Z2.ave(),M.ave(),T/2.0,iT,+1);
  
  auto tmp=
    multiFit<I+1>(corr,T,basePath,glbTMin,glbTMax,args...);
  
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
      
      grace_file_t plot(basePath/"glbFit.xmg");
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

enum OSTM{OS,TM};

string osTmTag(const OSTM& osTm)
{
  return ((osTm==OS)?"OS":"TM");
}

enum{TRANSVERSE,LONGITUDINAL};

template <bool TL>
struct ChannAnalysis
{
  const string chann;
  
  size_t tMinState0;
  
  size_t tMaxState0;
  
  size_t tMinState1;
  
  size_t tMaxState1;
  
  size_t tMinStates;
  
  size_t tMaxStates;
  
  size_t tSpecDensNorm;
  
  double EMinIntInGeV;
  
  double resWidthInGeV;
  
  double lambda;
  
  static constexpr size_t tMinBasis=1;
  
  size_t tMaxBasis;
  
  ChannAnalysis(const string& chann) :
    chann(chann)
  {
  }
  
  djack_t analyze(const djvec_t& fullCorr,
		  const vector<double>& fullCov,
		  const string& ensName,
		  const OSTM& osTm,
		  const size_t& T,
		  const djack_t& aInv,
		  const double& z) const
  {
    cout<<"*** "<<chann<<" ***"<<endl;
    
    const path baseOut=
      (path(::baseOut)/ensName/osTmTag(osTm)/chann);
    create_directories(baseOut);
    
    const Basis<Real> basis(T,tMinBasis,tMaxBasis,+1);
    const size_t nT=basis.nT;
    
    const auto fitPars=
      multiFit<0>(fullCorr,T,baseOut,tMinStates,tMaxStates,tMinState0,tMaxState0,"state0Fit.xmg",tMinState1,tMaxState1,"state1Fit.xmg");
    cout<<"state from multifit: "<<(fitPars[1].second*aInv).ave_err()<<endl;
    
    const double EMinInt=
      EMinIntInGeV/aInv.ave();
    const double EMaxInt=4;
    
    const djack_t normFactInGeV=
      sqr(z)*sqr(aInv)*aInv.ave()*24*sqr(M_PI)/cube(mTau);
    
    const double sigma=0;
    
    djvec_t corr(nT);
    
    MatrixX<Real> corrCov(nT,nT);
    for(size_t iT=0;iT<nT;iT++)
      {
	const size_t t=iT+tMinBasis;
	
	corr[iT]=fullCorr[t];
	for(size_t iS=0;iS<nT;iS++)
	  corrCov(iT,iS)=fullCov[t+(T/2+1)*(iS+tMinBasis)];
      }
    
    corr.ave_err().write(baseOut/"corrToReco.xmg");
    effective_mass(corr,T/2,0).ave_err().write(baseOut/"eff_mass_corrToReco.xmg");
    
    /// Preconditioner of the problem
    vector<Real> preco(nT);
    for(size_t iT=0;iT<basis.nT;iT++)
      preco[iT]=1/corr[iT].ave();
    
    const auto targetFunction=
      [sigma,amTau=mTau/aInv.ave()](const Real& E)
      {
	const Real x=E/amTau;
	
	return ((TL==TRANSVERSE)?
		Kt:Kl)(sigma,x);
      };
    
    const auto getEnlargedResonance=
      [aInv=aInv.ave()](const djack_t& Z2,
			const djack_t& M,
			const double& widthInGeV)
      {
	return [width=widthInGeV/aInv,
		Z2=Z2.ave(),
		M=M.ave()](const double& E)
	{
	  return Z2*gauss(M,width,E)/(2*M);
	};
      };
    
    
      // //parameters to fit
      // jack_fit_t fit;
      // const auto& [Z2Guess,MGuess]=fitPars[1];
      
      // djvec_t p(5);
      // p[0]=Z2Guess.ave();
      // p[1]=MGuess.ave();
      // p[2]=0.03;
      // p[3]=0.75*MGuess.ave();
      // p[4]=0.01;
      // cout<<"Guesses:"<<endl<<p.ave()<<endl;

      // fit.add_fit_par_limits(p[0],"Z2G",Z2Guess.ave_err(),Z2Guess.ave()/1.2,Z2Guess.ave()*1.2);
      // fit.add_fit_par_limits(p[1],"MG",MGuess.ave_err(),MGuess.ave()/1.2,MGuess.ave()*1.2);
      // fit.add_fit_par_limits(p[2],"WG",{p[2].ave(),0.001},0.01,0.3);
      // fit.add_fit_par_limits(p[3],"ThresC",{p[3].ave(),MGuess.err()},MGuess.ave()/3,3*MGuess.ave());
      // fit.add_fit_par_limits(p[4],"WC",{0.05,0.001},0,1);
      // const auto ttt=
      // 	[](const double& E,
      // 	   const double& t)
      // 	{
      // 	  return exp(-E *t)* (-E*E/t - (2 *E)/(t*t) - 2/(t*t*t));
      // 	};
      
      // const auto f=
      // 	[&basis,
      // 	 &EMaxInt,
      // 	 &EMinInt,
      // 	 &ttt](const auto& p,
      // 	       const double& t)
      // 	{
      // 	  return
      // 	    gslIntegrateFromTo([&basis,
      // 				t,
      // 				&p](const double& E)
      // 	    {
      // 	      const auto d=
      // 		p[0]*gauss(p[1],p[2],E)/(2*E);
      // 	      return d*basis(t-1,E);
      // 	    },std::min(EMinInt,p[1]-5*fabs(p[2])),std::max(EMaxInt,p[1]+5*fabs(p[2])))+
      // 	    p[4]*(ttt(4,t)-ttt(p[3],t))// gslIntegrateFromTo([&basis,
      // 	    // 			t](const double& E)
      // 	    // {
      // 	    //   const auto d=
      // 	    // 	E*E;
      // 	    //   return d*basis(t-1,E);
      // 	    // },p[3],EMaxInt);
      // 			       ;
      // 	};
      
      // for(size_t t=tMinStates;t<tMaxStates;t++)
      // 	fit.add_point(//numerical data
      // 		      [o=fullCorr[t]]
      // 		      (const vector<double>& p,
      // 		       const int& iel)
      // 		      {
      // 			return o[iel];
      // 		      },
      // 		      //ansatz
      // 		      [&f,t]
      // 		      (const vector<double>& fp,
      // 		       const int& iel)
      // 		      {
      // 			return f(fp,t);
      // 		      },
      // 		      //for covariance/error
      // 		      fullCorr[t].err());
      
      // // double initCh2=0;
      // // for(size_t iT=glbTMin;iT<glbTMax;iT++)
      // // 	{
      // // 	  const double a=in[iT].ave();
      // // 	  const double t=f(p.ave(),iT);
      // // 	  const double e=in[iT].err();
      // // 	  const double contr=sqr((a-t)/e);
	  
      // // 	  cout<<"contr=(sqr(("<<a<<"-"<<t<<")/"<<e<<")"<<endl;
      // // 	  initCh2+=contr;
      // // 	}
      // // cout<<"Initial ch2: "<<initCh2<<endl;
      
      // // fit_debug=true;
      // fit.fit();
      
      // // double finalCh2=0;
      // // for(size_t iT=glbTMin;iT<glbTMax;iT++)
      // // 	finalCh2+=sqr((in[iT].ave()-f(p.ave(),iT))/in[iT].err());
      // // cout<<"Final ch2: "<<finalCh2<<endl;

      // grace_file_t plot(baseOut+"glbFit2.xmg");
      // plot.write_vec_ave_err(effective_mass(fullCorr).ave_err());
      // plot.write_polygon([&](const double& x)
      // {
      // 	const auto u=
      // 	  [&p,&f](const double& x)
      // 	  {
      // 	    return jackCall(f,p,x);
      // 	  };
      // 	return effective_mass(u(x),u(x+1),x,T/2);
      // },tMinStates,tMaxStates,grace::BLUE);
      
      // cout<<p.ave_err()<<endl;
      // cout<<(p[2]*aInv).ave_err()<<endl;
    
    
    const auto& [Z2GroundState,MGroundState]=
      fitPars[1];
    const auto enlargedRes=
      getEnlargedResonance(Z2GroundState,MGroundState,resWidthInGeV);
      
    djvec_t cNorm(nT);
    for(size_t iT=0;iT<nT;iT++)
      {
	const double groundStateContr=
	  integrateUpToInfinity<double>([&basis,
					 iT,
					 &enlargedRes](const double& E)
	  {
	    return enlargedRes(E)*basis(iT,E);
	  });
	
	/// Correlator assuming that R(E)=E^N
	const Real cAss=
	  gslIntegrateFromTo([iT,&basis](const double& E)
	  {
	    return
	      get(
		  //exp(3.99*E)
		  pow(E,RE_exp)
		  
		  *basis(iT,E));
	  },MGroundState.ave()-resWidthInGeV/aInv.ave(),EMaxInt);
	
	cNorm[iT]=(fullCorr[iT]-groundStateContr)/get(cAss);
      }
    cNorm.ave_err().write(baseOut/"norm.xmg");
    
    // djvec_t cNorm2(nT);
    // for(size_t iT=0;iT<nT;iT++)
    //   {
    // 	const size_t t=iT+1;
    // 	const djack_t c=
    // 	  jackCall(f,p,t);
	
    // 	cNorm2[iT]=fullCorr[t]/c;
    //   }
    // cNorm2.ave_err().write(baseOut+"norm2.xmg");
    
    /// Time at which to normalize the correlator
    cout<<"tNorm: "<<tSpecDensNorm<<endl;
    
    /// Ansatz of the spectral density
    const auto specAns=
      [norm=cNorm[tSpecDensNorm].ave(),
       &enlargedRes// ,
       // p=p.ave()
       ](const Real& E)
      {
	return exp(3.99*E/2);
	// return p[0]*gauss(p[1],p[2],E)/(2*p[1])+
	//   p[4]*E*E*(E>p[3])*(E<4);
	return norm*pow(E,RE_exp)+enlargedRes(E);
      };
    
    grace_file_t specAnsPlot(baseOut/"specAns.xmg");
    specAnsPlot.write_line(specAns,EMinInt,EMaxInt);
    
    const TrambaccoFunctional tf=
      getTrambaccoFunctional<Real>(basis,
				   specAns,
				   targetFunction,
				   EMinInt,EMaxInt,corrCov,preco);
    
    /// Track the results
    djack_t R;
    VectorX<Real> g;
    
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
	
	Rprec*=normFactInGeV;
	
	jackCall([](double& R,
		    const Real& Rprec)
	{
	  R=get(Rprec);
	},R,Rprec);
      };
    
    // Performs the stability check
    grace_file_t stabilityPlot(baseOut/"stab.xmg");
    stabilityPlot.set_xaxis_logscale();
    double stabX;
    vector<pair<double,djack_t>> stab;
    for(double statOverSyst=64;statOverSyst>=1e-10;statOverSyst/=2)
      {
	reconstruct(statOverSyst);
	
	auto getA=
	  [&](const bool& withG)
	  {
	    return gslIntegrateFromTo([&targetFunction,
				       nT,
				       &g,
				       &basis,
				       &withG](const double& E)
	    {
	      Real s=0;
	      if(withG)
		for(size_t iT=0;iT<nT;iT++)
		  s+=g[iT]*basis(iT,E);
	      
	      return sqr(targetFunction(E)-s);
	    },EMinInt,EMaxInt);
	  };
	
	//A+B/statOverSyst;
	// const double lambda=1/(statOverSyst+1);
	const double Ag=getA(true);
	const double A0=getA(false);
	stabX=Ag/A0;
	stab.emplace_back(stabX,R);
      }
    
    // Final reconstruction
    reconstruct(exp(lambda));
    
    for(const auto& [x,y] : stab)
      stabilityPlot.write_ave_err(x,y.ave_err());
    stabilityPlot.write_constant_band(stabX-0.1,stabX+0.1,R);
    
    /// Print coefficients
    grace_file_t gPlot(baseOut/"g.xmg");
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
    grace_file_t recoPlot(baseOut/"reco.xmg");
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
    grace_file_t errPlot(baseOut/"err.xmg");
    errPlot.write_line([&](const Real& x)
    {
      const double reco=get(interpReco(x)*x);
      
      const Real EInGeV=x;
      const Real E=EInGeV/aInv.ave();
      
      const double targ=get(targetFunction(E)*x);
      
      return reco-targ;
    },EMinInt*aInv.ave(),EMaxInt*aInv.ave(),1000);
    
    /// Prints the reconstruction
    cout<<"Channel contribution: "<<smart_print(R)<<endl;
    
    return R;
  }
};

struct EnsDetails
{
  const double am;
  const size_t T;
  const size_t L;
  const string pref;
};

auto analyzeEns(const string& ensName,
		const OSTM& osTm,
		const string& sValue)
{
  const map<char,Real> zVlist{{'Z',0.697},{'B',0.706379},{'C',0.725404},{'D',0.744108},{'E',0.7663}};
  const map<char,Real> zAlist{{'Z',0.7357},{'B',0.74294},{'C',0.75830},{'D',0.77395},{'E',0.7921}};
  const map<char,ave_err_t> aList{{'B',{0.0795739,0.000132632}},{'C',{0.0682083,0.000134938}},{'D',{0.0569183,0.000115387}}};
  
  map<string,EnsDetails> ensDetails{{"B.72.64", {0.00072, 128, 64, "mix_fixed"}},
				    {"C.06.80", {0.00060, 160, 80, ""}},
				    {"C.06.112",{0.00060, 224, 112,"mix_fixed"}},
				    {"D.54.96", {0.00054, 192, 96, "mix"}}};
  
  const auto& eDet=ensDetails.at(ensName);
  // const double am=eDet.am;
  const size_t T=eDet.T;
  const size_t L=eDet.L;
  const string pref=eDet.pref;
  
  const double zOfV=
    (osTm==TM?zAlist:zVlist).at(ensName[0]);
  
  const double zOfA=
    (osTm==TM?zVlist:zAlist).at(ensName[0]);
  
  auto loadCorrCov=
    [&ensName,
     &T](const string& chann)
    {
      const string baseIn=
	"/home/francesco/QCD/LAVORI/GM3/";
      
      const auto [readT,corr,cov]=
	::loadCorrCov(baseIn+"/"+ensName+"/data/"+chann,false);
      
      if(T!=readT)
	CRASH("T=%zu, read %zu",T,readT);
      
      corr.ave_err().write(baseOut+".xmg");
      effective_mass(corr).ave_err().write(baseOut+"/"+ensName+"/"+chann+"eff_mass.xmg");
      
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
  
  // const djack_t aInv=
  //   estimateA(am,cP5P5_ll,T,baseOut+"/"+ensName+"/pion.xmg");
  const djack_t aInFm(gauss_filler_t{aList.at(ensName[0]),ensName[0]+1411325});
  const djack_t aInv=0.197/aInFm;
  cout<<"a^-1: "<<aInv.ave_err()<<" GeV^-1"<<endl;
  cout<<"a: "<<aInFm.ave_err()<<" fm"<<endl;
  
#define LOAD(TAG) \
  const auto [c ## TAG ## _ls,cov ## TAG ## _ls]=\
    loadCorrCov(pref + "_l_"+sValue+"_" + osTmTag(osTm) + "_" #TAG)
  
  LOAD(V0V0);
  LOAD(VKVK);
  LOAD(A0A0);
  LOAD(AKAK);
  LOAD(P5P5);
  
#undef LOAD
  
  const djack_t MPion=
    constant_fit(effective_mass(cP5P5_ll),30,60,baseOut+ensName+"/pion_fit.xmg");
  const djack_t MKaon=
    constant_fit(effective_mass(cP5P5_ls),30,60,baseOut+ensName+"/kaon_fit.xmg");
  
  const auto [Z2KaonStar,MKaonStar]=
    two_pts_fit(cVKVK_ls,T/2,22,35,baseOut+ensName+"/kaonStar_fit.xmg");
  const auto [Z2K1Star,MK1Star]=
    two_pts_fit(cAKAK_ls,T/2,20,24,baseOut+ensName+"/k1Star_fit.xmg");
  
  const double aPmin=2*M_PI/L;
  
  cout<<"MK: "<<smart_print(MKaon*aInv)<<" GeV"<<endl;
  cout<<"MK*: "<<smart_print(MKaonStar*aInv)<<" GeV"<<endl;
  cout<<"MK1*: "<<smart_print(MK1Star*aInv)<<" GeV"<<endl;
  
  const djack_t EKPi_nonint_FV=sqrt(sqr(MKaon)+sqr(aPmin))+sqrt(sqr(MPion)+sqr(aPmin));
  cout<<"E(K+Pi) non int FVE: "<<smart_print(EKPi_nonint_FV*aInv)<<" GeV"<<endl;
  
  auto readPars=
    [&ensName](auto& c,
	     const size_t n,
	     const size_t nMax)
    {
      raw_file_t input(ensName+"/pars.txt","r");
      
      for(size_t i=0;i<n;i++)
        input.read<string>();
      input.expect({c.chann});
      for(size_t i=n+1;i<nMax;i++)
        input.read<string>();
      
#define READ(TAG)					\
      input.expect(#TAG);				\
      for(size_t i=0;i<n;i++)				\
	input.read<string>();				\
      input.read(c.TAG);				\
      for(size_t i=n+1;i<nMax;i++)			\
	input.read<string>()				\
      
      READ(tMinState0);
      READ(tMaxState0);
      READ(tMinState1);
      READ(tMaxState1);
      READ(tMinStates);
      READ(tMaxStates);
      READ(EMinIntInGeV);
      READ(tMaxBasis);
      READ(tSpecDensNorm);
      READ(resWidthInGeV);
      READ(lambda);
      
#undef READ
    };
  
  ChannAnalysis<TRANSVERSE> VKVKchann("VKVK");
  ChannAnalysis<TRANSVERSE> AKAKchann("AKAK");
  ChannAnalysis<LONGITUDINAL> V0V0chann("V0V0");
  ChannAnalysis<LONGITUDINAL> A0A0chann("A0A0");
  readPars(VKVKchann,0,4);
  readPars(AKAKchann,1,4);
  readPars(V0V0chann,2,4);
  readPars(A0A0chann,3,4);
  
  const djack_t VKVKcontr=VKVKchann.analyze(cVKVK_ls,covVKVK_ls,ensName,osTm,T,aInv,zOfV);
  const djack_t AKAKcontr=AKAKchann.analyze(cAKAK_ls,covAKAK_ls,ensName,osTm,T,aInv,zOfA);
  const djack_t V0V0contr=V0V0chann.analyze(cV0V0_ls,covV0V0_ls,ensName,osTm,T,aInv,zOfV);
  const djack_t A0A0contr=A0A0chann.analyze(cA0A0_ls,covA0A0_ls,ensName,osTm,T,aInv,zOfA);
  
  const djack_t tot=VKVKcontr+AKAKcontr+V0V0contr+A0A0contr;
  
  cout<<"Total: "<<smart_print(tot)<<endl;
}

int main()
{
  set_njacks(50);
  if constexpr(std::is_same_v<PrecFloat,Real>)
    PrecFloat::setDefaultPrecision(128);
  
  for(const string ensName : {"B.72.64","C.06.80","C.06.112","D.54.96"})
    for(const OSTM& osTm : {OS,TM})
      for(const string sValue : {"s1"})
	{
	  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
	  cout<<"////////////////////////////// "<<ensName<<" "<<osTmTag(osTm)<<" ///////////////////////////////"<<endl;
	  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
	  
	  analyzeEns(ensName,osTm,sValue);
	}
  
  return 0;
}
