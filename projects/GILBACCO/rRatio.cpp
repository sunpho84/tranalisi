#include <tranalisi.hpp>

#include <loadCorrCov.hpp>

#include <trambacco.hpp>

using Real=double;
using Jack=jack_t<Real>;
using Jvec=vmeas_t<Jack>;
using Vect=Vector<Real,Dynamic>;
using Matr=Matrix<Real,Dynamic,Dynamic>;
constexpr Real tol=1e-8;

/// Long distance behaviour of spectral density: E^N
const Real RE_exp=2;

/// Minimal time from which to consider the correlator
const size_t tMin=1;

/// Lattice extention
size_t T;

/// Width of the distribution at EMin
const Real sigmaInGeVEMin=0.200;

/// Scale sigma with energy
constexpr bool scaleSigma=true;

/// Physical pion mass
constexpr Real MPiPhys=0.1350;

/// Physical pion decay constant
constexpr Real fPiPhys=0.1304;

/// Base directory to read data
const string baseIn="/home/francesco/QCD/LAVORI/GM3";

/// Base directory to write data
const string baseOut="/home/francesco/QCD/LAVORI/TRAMBACCO";

/// Loads a correlator avoiding to compute covariance, dropping returned dummy vector
djvec_t loadCorr(const string& ensName,
		 const string& chann)
{
  const auto [T,corr,cov]=
    loadCorrCov(baseIn+"/"+ensName+"/data/"+chann,false);
  
  return corr;
}

/// Analyse a given ensemble
void analyzeEns(grace_file_t& glbRplot,
		const string& ensName)
{
  // Don't move, needed to set T
  const djvec_t cP5P5=
    loadCorr(ensName,"ll_TM_P5P5");
  cP5P5.ave_err().write(baseOut+"/"+ensName+"/P5P5.xmg");
  
  const map<char,Real> zVlist{{'Z',0.697},{'B',0.706379},{'C',0.725404},{'D',0.744108},{'E',0.7663}};
  const map<char,Real> zAlist{{'Z',0.7357},{'B',0.74294},{'C',0.75830},{'D',0.77395},{'E',0.7921}};
  const map<char,Real> aMlist{{'Z',0.00077},{'B',0.00072},{'C',0.0006},{'D',0.00054},{'E',0.00044}};
  const map<char,size_t> tMaxList{{'Z',26},{'B',30},{'C',36},{'D',40},{'E',35}};
  const Real z=zAlist.at(ensName[0]);
  const Real am=aMlist.at(ensName[0]);
  // const size_t tMax=T/2+1;
  const size_t tMax=tMaxList.at(ensName[0]);
  const size_t nT=tMax-tMin;
  
  // Lattice spacing determination
  djack_t aMPi,Z2Pi;
  two_pts_fit(Z2Pi,aMPi,cP5P5,T/2,30,40,baseOut+"/"+ensName+"/pion.xmg");
  const djack_t aFPi=sqrt(Z2Pi)*2*am/(aMPi*sinh(aMPi));
  cout<<"aMPi: "<<aMPi.ave_err()<<endl;
  cout<<"aFPi: "<<aFPi.ave_err()<<endl;
  const djack_t xi=aMPi/aFPi;
  cout<<"Xi: "<<xi.ave_err()<<" phys: "<<MPiPhys/fPiPhys<<endl;
  cout<<"Z: "<<z<<endl;
  
  const djvec_t cA0A0=
    loadCorr(ensName,"ll_TM_A0A0");
  cA0A0.ave_err().write(baseOut+"/"+ensName+"/A0A0.xmg");
  const djack_t aMPiFromA0A0=constant_fit(effective_mass(cA0A0),8,T/2-3,baseOut+"/"+ensName+"/effA0A0.xmg");
  
  cout<<aMPi.ave_err().significativity()<<" "<<aMPiFromA0A0.ave_err().significativity()<<endl;

  // {
  // const djvec_t cV0P5=
  //   loadCorr(ensName,"ll_TM_V0P5");
  // cV0P5.ave_err().write(baseOut+"/"+ensName+"/corr_V0P5.xmg");
  
  // const djvec_t cV0P5der=
  //   symmetric_derivative(cV0P5);
  
  // cV0P5der.ave_err().write(baseOut+"/"+ensName+"/corr_V0P5_der.xmg");
  
  // const djvec_t mPCACEff=
  //   cV0P5der/(2*cP5P5);
  
  // const djack_t mPCAC=
  //   constant_fit(mPCACEff,14,T/2-1,baseOut+"/"+ensName+"/eff_mPCAC.xmg");
  
  // cout<<"mPCAC: "<<mPCAC.ave_err()<<endl;
  // }
  
  /// Inverse lattice spacing, in GeV
  const djack_t aInv=
    fPiPhys/aFPi;
  cout<<"a^-1: "<<smart_print(aInv)<<" GeV"<<endl;
  
  /// Lattice spacing in fm
  const djack_t aInFm=
    0.197/aInv;
  cout<<ensName<<" a: "<<smart_print(aInFm)<<" fm"<<endl;
  
  /// Loads the correlation function of which to compute R
  const auto [readT,cVKVK,covVKVK]=
    loadCorrCov(ensName,"ll_TM_VKVK");
  
  if(T==0)
    T=readT;
  else
    if(T!=readT)
      CRASH("read size %zu different from stored size %zu",readT,T);
  
  cVKVK.ave_err().write(baseOut+"/"+ensName+"/VKVK.xmg");
  effective_mass(cVKVK).ave_err().write(baseOut+"/"+ensName+"/rho.xmg");
  
  /// Prepare the renormalized correlator
  Jvec corr(nT);
  vector<Real> corrCov(nT*nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      const Real fact=sqr(z);//*12*sqr(M_PI)*4/9
      corr[iT]=cVKVK[iT+tMin]*fact;
      for(size_t iS=0;iS<nT;iS++)
	corrCov[iT+nT*iS]=covVKVK[(iT+tMin)+(T/2+1)*(iS+tMin)]*sqr(fact);
    }
  corr.ave_err().write(baseOut+"/"+ensName+"/rhoPreco.xmg");
  
  /// Preconditioner of the problem
  vector<Real> preco(nT);
  for(size_t iT=0;iT<nT;iT++)
    preco[iT]=corr[iT].ave();
  
  grace_file_t RPlot(baseOut+"/"+ensName+"/R.xmg");
  const Real EMinInGeV=0.270,EMaxInGeV=4;
  const Real EMin=EMinInGeV/aInv.ave();
  const Real EMax=EMaxInGeV/aInv.ave();
  const Real EMaxInt=EMax*4;
  
  /// Normalization of the correlation function, when assuming that R(E)=E^N
  Jvec norm(nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      /// Correlator assuming that R(E)=E^N
      const double cAss=
	gslIntegrateFromTo([iT](const Real& E)
	{
	  return pow(E,RE_exp)*b(iT,E);
	},EMin,EMaxInt,tol);
      
      norm[iT]=corr[iT]/cAss;
    }
  norm.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  /// Time at which to normalize the correlator: we choose 2 GeV^-1
  const size_t tNorm=round(2*aInv.ave());
  cout<<"tNorm: "<<tNorm<<endl;
  
  /// Plots the value of sigma as a function of the energy
  grace_file_t sigmaPlot(baseOut+"/"+ensName+"/sigma.xmg");
  
  struct Res
  {
    Jack R;
    double sigma;
    double syst;
    double stat;
  };
  
  vector<tuple<Real,Jack>> RofE;
  vector<double> sigmaOfE;
  vector<double> systOfE,statOfE;
  
  /// Holds the reconstruction of the target function
  vector<function<double(double)>> recoOfE;
  
  for(Real sigma,EStar=EMin;EStar<EMax;EStar+=sigma/10)
    {
      /// Computes Estar in GeV
      const Real EStarInGeV=EStar*aInv.ave();
      sigma=sigmaInGeVEMin/aInv.ave()*(scaleSigma?sqrt(EStar/EMin):1);
      
      /// Compute sigma in GeV
      const Real sigmaInGeV=sigma*aInv.ave();
      sigmaOfE.emplace_back(sigmaInGeV);
      
      cout<<"E: "<<EStarInGeV<<" GeV, Sigma: "<<sigmaInGeV<<" GeV"<<endl;
      
      sigmaPlot.write_xy(EStarInGeV,sigmaInGeV);
      
      /// Target function
      const auto targ=
	[EStar,sigma](const Real& E)
	{
	  return gauss(EStar,sigma,E)/sqr(E);
	};
      
      Vect e(nT);
      Matr f(nT,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  e(iT)=
	    gslIntegrateFromTo([targ,iT](const double& E)
	    {
	      return targ(E)*pow(E,2*RE_exp)*b(iT,E);
	    },EMin,EMaxInt,tol);
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      gslIntegrateFromTo([iT,iS](const double& E)
	      {
		return pow(E,2*RE_exp)*b(iT,E)*b(iS,E);
	      },EMin,EMaxInt,tol);
	}
      
      
      
      // for(size_t iT=0;iT<nT;iT++)
      // 	{
      // 	  e(iT)=
      // 	    gslIntegrateUpToInfinity([targ,iT](const double& E)
      // 	    {
      // 	      return targ(E)*pow(E,2*N)*b(iT,E);
      // 	    },EMin,tol);
	  
      // 	  for(size_t iS=iT;iS<nT;iS++)
      // 	    f(iT,iS)=f(iS,iT)=
      // 	      gslIntegrateUpToInfinity([iT,iS](const double& E)
      // 	      {
      // 		return pow(E,2*N)*b(iT,E)*b(iS,E);
      // 	      },EMin,tol);
      // 	}
      
      const Real c2=
	sqr(norm[tNorm].ave());
      
      VectorXd L(nT);
      MatrixXd QStat(nT,nT),QSyst(nT,nT);
      
      /// Prepares the functional components
      for(size_t iT=0;iT<nT;iT++)
	{
	  L(iT)=c2*e(iT)/preco[iT];
	  for(size_t iS=0;iS<nT;iS++)
	    {
	      constexpr double corrFading=1;
	      const size_t dist=abs((int)iS-(int)iT);
	      QSyst(iT,iS)=c2*f(iT,iS)/(preco[iT]*preco[iS]);
	      QStat(iT,iS)=corrCov[iT+nT*iS]*pow(corrFading,dist)/(preco[iT]*preco[iS]);
	    }
	}
      
      /// Track the results
      Vect g;
      Jack R;
      
      const auto reconstruct=
	[&](const Real& statOverSyst=1)
	{
	  MatrixXd Q=
	    QSyst+QStat/statOverSyst;
	  
	  g=Q.inverse()*L;
	  
	  R=0;
	  for(size_t iT=0;iT<nT;iT++)
	    R+=g[iT]*corr[iT]/preco[iT];
	};
      
      // Performs the stability check
      grace_file_t stabilityPlot(baseOut+"/"+ensName+"/stab"+to_string(EStarInGeV)+".xmg");
      vector<pair<double,djack_t>> stab;
      for(Real statOverSyst=512;statOverSyst>=1e-4;statOverSyst/=2)
	{
	  reconstruct(statOverSyst);
	  stab.emplace_back(log(statOverSyst),R);
	}
      
      // Final reconstruction
      const double lambdaExp=-2.5;
      reconstruct(exp(lambdaExp));
      
      for(const auto& [x,y] : stab)
	stabilityPlot.write_ave_err(x,y.ave_err());
      stabilityPlot.write_constant_band(lambdaExp-0.1,lambdaExp+0.1,R);
      
      systOfE.push_back(sqrt(g.transpose()*QSyst*g));
      statOfE.push_back(sqrt(g.transpose()*QStat*g));
      
      /// Print coefficients
      grace_file_t gPlot(baseOut+"/"+ensName+"/g"+to_string(EStarInGeV)+".xmg");
      for(size_t iT=0;iT<nT;iT++)
	gPlot.write_xy(iT,g[iT]);
      
      /// Estimates the cancellation
      // Jack maxCoeff;
      // maxCoeff=0;
      // {
      // 	Jack r;
      // 	r=0;
      // 	for(size_t iT=0;iT<nT;iT++)
      // 	  {
      // 	    const Jack contr=g[iT]*corr[iT]/preco[iT];
      // 	    r+=contr;
      // 	    maxCoeff=max(abs(contr),maxCoeff);
      // 	  }
      // 	const Jack cancel=maxCoeff/r;
      // 	cout<<EStar<<" "<<cancel.ave_err()<<endl;
      // }
      // Real errEst=0;
      // for(size_t iT=0;iT<nT;iT++)
      // 	for(size_t iS=0;iS<nT;iS++)
      // 	  errEst+=gRes[iT].ave()*gRes[iS].ave()*v[iS+nT*iT]*(iS==iT || 1)/(preco[iT]*preco[iS]);
      // cout.precision(16);
      // cout<<R.err()<<" "<<sqrt(errEst)<<endl;
      
      const auto interpReco=
	[g,&aInv,&nT,&preco,&trambacco](const double& EInGeV)
      {
	const Real E=EInGeV/aInv.ave();
	Real s;
	s=0;
	for(size_t iT=0;iT<nT;iT++)
	  {
	    const Real k=trambacco.b(iT,E)/preco[iT];
	    s+=g[iT]*k;
	  }
	s*=sqr(E);
	
	return s;
      };
      
      recoOfE.emplace_back(interpReco);
      
      // Plots the reconstruction of the target function
      grace_file_t recoPlot(baseOut+"/"+ensName+"/reco"+to_string(EStarInGeV)+".xmg");
      recoPlot.write_line(interpReco,EMin*aInv.ave(),EMax*aInv.ave());
      
      recoPlot.write_line([&](const Real& EInGeV)
      {
	const Real E=EInGeV/aInv.ave();
	return targ(E)*sqr(E);
      },EMin*aInv.ave(),EMax*aInv.ave());
      
      RofE.emplace_back(EStarInGeV,R);
    }
  
  glbRplot.set_legend(ensName);
  
  const auto interpR=
    [&RofE](const Real& E)
  {
    vector<pair<Real,size_t>> distPos(RofE.size());
    for(size_t iE=0;iE<RofE.size();iE++)
      distPos[iE]={fabs(get<0>(RofE[iE])-E),iE};
    sort(distPos.begin(),distPos.end());
    
    vector<double> x(3);
    Jvec y(3);
    for(size_t iE=0;iE<3;iE++)
      tie(x[iE],y[iE])=RofE[distPos[iE].second];
    
    const Jvec pars=poly_fit(x,y,2);
    
    return poly_eval(pars,E);
  };
  
  for(grace_file_t* plot : {&RPlot,&glbRplot})
      plot->write_polygon(interpR,get<0>(RofE.front()),get<0>(RofE.back()));
  
  grace_file_t biConv(baseOut+"/"+ensName+"/biConv.xmg");
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
  //     grace_file_t("").write_polygon([&](const double& E)
  //     {
  // 	return abs(std::get<1>(recoOfE[iE])(E)*sqr(E)*interpR(E));
  // }, double xmin, double xmax, grace::color_t col)
      
      const double x=std::get<0>(RofE[iE]);
      
      grace_file_t("/tmp/"+to_string(x)+".xmg").write_polygon([&](const double& E)
      {
	return (recoOfE[iE](E)/aInv.ave()-gauss(x,sigmaOfE[iE],E))*interpR(E);
      },EMinInGeV,EMaxInGeV);
      
      const double syst=
	gslIntegrateUpToInfinity([&](const double& E)
      {
	return sqr((recoOfE[iE](E)/aInv.ave()-gauss(x,sigmaOfE[iE],E))*interpR(E).ave());
      },EMinInGeV,1e-2);
      biConv.write_xy(x,sqrt(syst));
    }
  biConv.set_all_colors(grace::BLACK);
  
  biConv.new_data_set();
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
      const double x=std::get<0>(RofE[iE]);
      const double smooth=
	gslIntegrateUpToInfinity([&](const double& E)
      {
	return gauss(x,sigmaOfE[iE],E)*interpR(E).ave();
      },EMinInGeV,1e-2);
      biConv.write_xy(x,smooth);
    }
  biConv.set_all_colors(grace::MAGENTA);
  
  biConv.new_data_set();
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
      const double x=std::get<0>(RofE[iE]);
      biConv.write_xy(x,interpR(x).err());
    }
  biConv.set_all_colors(grace::RED);
  
  biConv.new_data_set();
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
      const double x=std::get<0>(RofE[iE]);
      biConv.write_xy(x,systOfE[iE]);
    }
  biConv.set_all_colors(grace::BLUE);
  
  biConv.new_data_set();
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
      const double x=std::get<0>(RofE[iE]);
      biConv.write_xy(x,statOfE[iE]);
    }
  biConv.set_all_colors(grace::GREEN4);
  
  biConv.new_data_set();
  for(size_t iE=0;iE<recoOfE.size();iE++)
    {
      const auto [x,y]=RofE[iE];
      biConv.write_xy(x,y.ave());
    }
  biConv.set_all_colors(grace::VIOLET);
  
  
  RPlot.write_line([c=norm[tNorm].ave()](const Real& x)
  {
    return c;
  },EMin*aInv.ave(),EMax*aInv.ave(),grace::BLACK);
}

vector<string> readEnsList()
{
  vector<string> ensList;
  ifstream ensListFile("ensList.txt");
  string s;
  while(ensListFile>>s)
    ensList.push_back(s);
  
  return ensList;
}

int main()
{
  set_njacks(50);
  
  grace_file_t RPlot(baseOut+"/R.xmg");
  RPlot.set_color_scheme({grace::BLACK,grace::RED,grace::GREEN4,grace::BLUE,grace::VIOLET});
  
  const vector<string> ensList=readEnsList();
  
  // for(const char* ensName : {"B.72.64","B.72.96","C.06.80","C.06.112","D.54.96"})
  for(const string& ensName : ensList)
  // // for(const char* ensName : {"C.06.80"})
  // for(const char* ensName : {"B.72.64"})
    analyzeEns(RPlot,ensName);
  
  return 0;
}
