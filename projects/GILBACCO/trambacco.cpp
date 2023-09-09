#include <tranalisi.hpp>

using Real=double;
using Jack=jack_t<Real>;
using Jvec=vmeas_t<Jack>;
using Vect=Vector<Real,Dynamic>;
using Matr=Matrix<Real,Dynamic,Dynamic>;
constexpr Real tol=1e-8;

/// Long distance behaviour of spectral density: E^N
const Real N=2;

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

/// Loads a correlator and compute covariance
///
/// We don't keep in count the correlation of the data, so the
/// covariance matrix as generated from the boostrap is
/// underestimated. So we pass to the correlation matrix, and then we
/// renormalize with the diagonal part taken from jackknive estimate
/// of the error
auto loadCorrCov(const string& ensName,
		 const string& chann,
		 const bool computeCov=true)
{
  /// File needed to read data
  raw_file_t fin(baseIn+"/"+ensName+"/data/"+chann,"r");
  
  /// Gets the number of confs
  const size_t nConfs=fin.bin_read<size_t>();
  
  /// Read and drops the number of hits
  [[maybe_unused]] const size_t nHits=
		     fin.bin_read<size_t>();
  
  // Read the time extention
  T=fin.bin_read<size_t>();
  
  /// Index to access the raw data
  index_t iRaw({{"nConfs",nConfs},{"T/2+1",T/2+1}});
  
  /// Raw data
  vector<double> raw(iRaw.max());
  fin.bin_read(raw);
  
  /// Performs the jackknife resampling of the correlator
  djvec_t res(T/2+1);
  jackknivesFill(nConfs,
		 [&](const size_t& iConf,
		     const size_t& iJack,
		     const double& w)
		 {
		   for(size_t t=0;t<=T/2;t++)
		     res[t][iJack]+=w*raw[iRaw({iConf,t})];
		 });
  res.clusterize((double)nConfs/njacks);
  
  /// Covariance matrix
  vector<double> covMatr((T/2+1)*(T/2+1),0.0);
  
  if(computeCov)
    {
      /// Random generator to create the boostrap sapmple
      gen_t gen(232452);
      
      /// Number of boostrap to be used
      const size_t nBoots=1000;
      
      /// Bootstrap average of the correlator
      vector<double> ave(T/2+1,0.0);
      
      /// Bootstrap sample
      vector<double> bootSamp(T/2+1,0.0);
      
      for(size_t iBoot=0;iBoot<nBoots;iBoot++)
	{
	  /// Reset the bootstrap
	  for(auto& bi : bootSamp)
	    bi=0.0;
	  
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    {
	      /// Draw a configuration
	      const size_t jConf=
		gen.get_int(0,nConfs);
	      
	      // Fills the bootstrap
	      for(size_t t=0;t<=T/2;t++)
		bootSamp[t]+=raw[iRaw({jConf,t})];
	    }
	  
	  /// Normalize the bootstrap
	  for(size_t t=0;t<=T/2;t++)
	    bootSamp[t]/=nConfs;
	  
	  /// Increment the average and compute correlation
	  for(size_t t=0;t<=T/2;t++)
	    {
	      ave[t]+=bootSamp[t];
	      for(size_t s=0;s<=T/2;s++)
		covMatr[s+(T/2+1)*t]+=bootSamp[t]*bootSamp[s];
	    }
	}
      
      // Normalize average and covariance
      for(size_t t=0;t<=T/2;t++)
	{
	  ave[t]/=nBoots;
	  for(size_t s=0;s<=T/2;s++)
	    covMatr[s+(T/2+1)*t]/=nBoots;
	}
      
      // Subtract product of averages from covariance
      for(size_t t=0;t<=T/2;t++)
	for(size_t s=0;s<=T/2;s++)
	  covMatr[t+(T/2+1)*s]-=ave[s]*ave[t];
      
      // Pass to the correlation matrix
      for(size_t t=0;t<=T/2;t++)
	{
	  for(size_t s=t+1;s<=T/2;s++)
	    {
	      covMatr[t+(T/2+1)*s]/=sqrt(covMatr[t+(T/2+1)*t])*sqrt(covMatr[s+(T/2+1)*s]);
	      covMatr[s+(T/2+1)*t]=covMatr[t+(T/2+1)*s];
	    }
	  covMatr[t+(T/2+1)*t]=1;
	}
      
      // Normalize back to the covariance matrix
      for(size_t t=0;t<=T/2;t++)
	{
	  const double d=res[t].err();
	  
	  for(size_t s=0;s<=T/2;s++)
	    {
	      covMatr[t+(T/2+1)*s]*=d;
	      covMatr[s+(T/2+1)*t]*=d;
	    }
	}
    }
  
  return make_tuple(res,covMatr);
}

/// Loads a correlator avoiding to compute covariance, dropping returned dummy vector
djvec_t loadCorr(const string& ensName,
		 const string& chann)
{
  const auto [corr,cov]=
    loadCorrCov(ensName,chann,false);
  
  return corr;
}

/// Factorial of n
constexpr int fact(const int& n)
{
  int r=1;
  
  for(int i=2;i<=n;i++)
    r*=i;
  
  return r;
}

/// Gaussian centered in x0 with standard deviation sigma, normalized to 1
Real gauss(const Real& x0,
	   const Real& sigma,
	   const Real& x)
{
  return exp(-sqr((x-x0)/sigma)/2)/(sqrt(2*M_PI)*sigma);
}

/// Basis function
///
/// The value of "it" is the element of the function, so we need to
/// increase by tmin to get the time
Real b(const size_t& it,
       const Real& E)
{
  /// Refered time
  const Real t=
    it+tMin;
  
  return exp(-t*E)+exp(-(T-t)*E);
}

/// Analyse a given ensemble
void analyzeEns(grace_file_t& glbRplot,
		const string& ensName)
{
  // Don't move, needed to set T
  const djvec_t cP5P5=
    loadCorr(ensName,"ll_TM_P5P5");
  
  const map<char,Real> zVlist{{'B',0.706379},{'C',0.725404},{'D',0.744108},{'E',0.76}};
  const map<char,Real> zAlist{{'B',0.74294},{'C',0.75830},{'D',0.77395},{'E',0.79}};
  const map<char,Real> aMlist{{'B',0.00072},{'C',0.0006},{'D',0.00054},{'E',0.00044}};
  const map<char,size_t> tMaxList{{'B',30},{'C',36},{'D',40},{'E',35}};
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
  
  /// Inverse lattice spacing, in GeV
  const djack_t aInv=
    fPiPhys/aFPi;
  cout<<"a^-1: "<<smart_print(aInv)<<" GeV"<<endl;
  
  /// Lattice spacing in fm
  const djack_t aInFm=
    0.197/aInv;
  cout<<ensName<<" a: "<<smart_print(aInFm)<<" fm"<<endl;
  
  /// Loads the correlation function of which to compute R
  const auto [cVKVK,covVKVK]=
    loadCorrCov(ensName,"ll_TM_VKVK");
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
  const Real EMin=EMinInGeV/aInv.ave(),EMax=EMaxInGeV/aInv.ave();
  
  /// Normalization of the correlation function, when assuming that R(E)=E^N
  Jvec norm(nT);
  for(size_t iT=0;iT<nT;iT++)
    {
      /// Correlator assuming that R(E)=E^N
      const double cAss=
	gslIntegrateUpToInfinity([iT](const Real& E)
	{
	  return pow(E,N)*b(iT,E);
	},EMin,tol);
      
      norm[iT]=corr[iT]/cAss;
    }
  norm.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  /// Time at which to normalize the correlator: we choose 2 GeV^-1
  const size_t tNorm=round(2*aInv.ave());
  cout<<"tNorm: "<<tNorm<<endl;
  
  /// Plots the value of sigma as a function of the energy
  grace_file_t sigmaPlot(baseOut+"/"+ensName+"/sigma.xmg");
  
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
	    gslIntegrateUpToInfinity([targ,iT](const double& E)
	    {
	      return targ(E)*pow(E,2*N)*b(iT,E);
	    },EMin,tol);
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      gslIntegrateUpToInfinity([iT,iS](const double& E)
	      {
		return pow(E,2*N)*b(iT,E)*b(iS,E);
	      },EMin,tol);
	}
      
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
      
      auto reconstruct=
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
      const double lambdaExp=-3;
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
	[g,&aInv,&nT,&preco](const double& EInGeV)
      {
	const Real E=EInGeV/aInv.ave();
	Real s;
	s=0;
	for(size_t iT=0;iT<nT;iT++)
	  {
	    const Real k=b(iT,E)/preco[iT];
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

int main()
{
  
  set_njacks(50);
  
  grace_file_t RPlot(baseOut+"/R.xmg");
  
  for(const char* ensName : {"B.72.64","B.72.96","C.06.80","C.06.112","D.54.96"})
  // // for(const char* ensName : {"C.06.80"})
  // for(const char* ensName : {"B.72.64"})
    analyzeEns(RPlot,ensName);
  
  return 0;
}
