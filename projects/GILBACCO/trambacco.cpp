#include <tranalisi.hpp>

const double N=2;
const size_t tMin=1;
size_t T;
const size_t tNorm=5;
const double sigmaInGeVEMin=0.200;
constexpr bool scaleSigma=true;
constexpr double corrFading=0.9;
constexpr double MPiPhys=0.1350,fPiPhys=0.1304;

const string baseIn="/home/francesco/QCD/LAVORI/GM3";
const string baseOut="/home/francesco/QCD/LAVORI/TRAMBACCO";

djvec_t loadCorr(const string& ensName,
		 const string& chann)
{
  raw_file_t fin(baseIn+"/"+ensName+"/data/"+chann,"r");
  const size_t nConfs=fin.bin_read<size_t>();
  [[maybe_unused]] const size_t nHits=fin.bin_read<size_t>();
  T=fin.bin_read<size_t>();
  
  index_t iRaw({{"nConfs",nConfs},{"T/2+1",T/2+1}});
  vector<double> raw(iRaw.max());
  fin.bin_read(raw);
  
  djvec_t corr(T/2+1);
  jackknivesFill(nConfs,
		 [&](const size_t& iConf,
		     const size_t& iJack,
		     const double& w)
		 {
		   for(size_t t=0;t<=T/2;t++)
		     corr[t][iJack]+=w*raw[iRaw({iConf,t})];
		 });
  corr.clusterize((double)nConfs/njacks);
  
  return corr;
}

constexpr int fact(const int& n)
{
  int r=1;
  
  for(int i=2;i<=n;i++)
    r*=i;
  
  return r;
}

double gauss(const double& x0,
	     const double& sigma,
	     const double& x)
{
  return exp(-sqr((x-x0)/sigma)/2)/(sqrt(2*M_PI)*sigma);
}

double b(const size_t& it,
	 const double& E)
{
  const double t=it+tMin;
  
  return exp(-t*E)+exp(-(T-t)*E);
}

void analyzeEns(grace_file_t& glbRplot,
		const string& ensName)
{
  const map<char,double> zVlist{{'B',0.706379},{'C',0.725404},{'D',0.744108},{'E',0.76}};
  const map<char,double> zAlist{{'B',0.74294},{'C',0.75830},{'D',0.77395},{'E',0.79}};
  const map<char,double> aMlist{{'B',0.00072},{'C',0.0006},{'D',0.00054},{'E',0.00044}};
  const map<char,size_t> tMaxList{{'B',30},{'C',36},{'D',40},{'E',35}};
  const double z=zAlist.at(ensName[0]);
  const double am=aMlist.at(ensName[0]);
  const size_t tMax=tMaxList.at(ensName[0]);
  const size_t nT=tMax-tMin;
  
  const djvec_t cP5P5=
    loadCorr(ensName,"ll_TM_P5P5");
  
  // Lattice spacing study
  djack_t aMPi,Z2Pi;
  two_pts_fit(Z2Pi,aMPi,cP5P5,T/2,30,40,baseOut+"/"+ensName+"/pion.xmg");
  const djack_t aFPi=sqrt(Z2Pi)*2*am/(aMPi*sinh(aMPi));
  cout<<"aMPi: "<<aMPi.ave_err()<<endl;
  cout<<"aFPi: "<<aFPi.ave_err()<<endl;
  const djack_t xi=aMPi/aFPi;
  cout<<"Xi: "<<xi.ave_err()<<" phys: "<<MPiPhys/fPiPhys<<endl;
  
  const djack_t aInv=
    fPiPhys/aFPi;
  const djack_t aInFm=
    0.197/aInv;
  cout<<"a^-1: "<<smart_print(aInv)<<" GeV"<<endl;
  cout<<ensName<<" a: "<<smart_print(aInFm)<<" fm"<<endl;
  
  const djvec_t cVKVK=
    loadCorr(ensName,"ll_TM_VKVK");
  effective_mass(cVKVK).ave_err().write(baseOut+"/"+ensName+"/rho.xmg");
  djvec_t corr(nT);
  for(size_t iT=0;iT<nT;iT++)
    corr[iT]=cVKVK[iT+tMin]*sqr(z);
  
  vector<double> preco(nT);
  for(size_t iT=0;iT<nT;iT++)
    preco[iT]=corr[iT].ave();
  
  // grace_file_t plot(baseOut+"/"+ensName+"/p.xmg");
  // plot.write_line([](const double& x){return ;},-10,10);
  
  grace_file_t RPlot(baseOut+"/"+ensName+"/R.xmg");
  const double EMin=2*aMPi.ave(),EMax=4/aInv.ave();
  
  djvec_t norm(nT);
  for(size_t iT=0;iT<nT;iT++)
    norm[iT]=corr[iT]/gslIntegrateFromTo([t=iT+tMin](const double& E)
    {
      return pow(E,N)*exp(-E*t);
    },EMin,EMax);
  
  norm.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  grace_file_t sigmaPlot(baseOut+"/"+ensName+"/sigma.xmg");
  
  double EStar=EMin;
  vector<tuple<double,djack_t>> RofE;
  while(EStar<EMax)
    {
      const double EStarInGeV=EStar*aInv.ave();
      const double sigma=sigmaInGeVEMin/aInv.ave()*(scaleSigma?sqrt(EStar/EMin):1);
      const double sigmaInGeV=sigma*aInv.ave();
      cout<<"E: "<<EStarInGeV<<" GeV, Sigma: "<<sigmaInGeV<<" GeV"<<endl;
      
      sigmaPlot.write_xy(EStarInGeV,sigmaInGeV);
      
      /// Target function
      const auto targ=
	[EStar,sigma](const double& E)
	{
	  return gauss(EStar,sigma,E)/sqr(E);
	};
      
      // const double D=
      // 	gslIntegrateUpToInfinity([targ](const double& E)
      // 	{
      // 	  return sqr(targ(E))*pow(E,2*N);
      // 	},EMin);
      
      VectorXd e(nT);
      MatrixXd f(nT,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  e(iT)=
	    gslIntegrateUpToInfinity([targ,iT](const double& E)
	    {
	      return targ(E)*pow(E,2*N)*b(iT,E);
	    },EMin);
	  // cout<<iT<<" "<<e(iT)<<endl;
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      gslIntegrateUpToInfinity([iT,iS](const double& E)
	      {
		return pow(E,2*N)*b(iT,E)*b(iS,E);
	      },EMin);
	}
      
      /// Computes the covariance matrix
      MatrixXd v(nT,nT);
      for(size_t iT=0;iT<nT;iT++)
	for(size_t iS=iT;iS<nT;iS++)
	  v(iT,iS)=v(iS,iT)=cov(corr[iT],corr[iS]);
      
      const double c2=
	sqr(norm[tNorm].ave());
      
      VectorXd L(nT);
      MatrixXd QStat(nT,nT),QSyst(nT,nT);
      
      /// Prepares the functional components
      for(size_t iT=0;iT<nT;iT++)
	{
	  L(iT)=c2*e(iT)/preco[iT];
	  for(size_t iS=0;iS<nT;iS++)
	    {
	      const size_t dist=abs((int)iS-(int)iT);
	      QSyst(iT,iS)=c2*f(iT,iS)/(preco[iT]*preco[iS]);
	      QStat(iT,iS)=v(iS,iT)*pow(corrFading,dist)/(preco[iT]*preco[iS]);
	    }
	}
      
      /// Track the results
      VectorXd g;
      djack_t R;
      
      auto reconstruct=
	[&](const double& statOverSyst=1)
	{
	  MatrixXd Q=
	    QSyst+QStat/statOverSyst;
	  
	  g=Q.inverse()*L;
	  
	  R=0;
	  for(size_t iT=0;iT<nT;iT++)
	    R+=g[iT]*corr[iT]/preco[iT];
	  
	  // const double syst=g.transpose()*QSyst*g;
	  // const double stat=(double)(g.transpose()*QStat*g)/statOverSyst;
	  
	  // cout<<statOverSyst<<" "<<sqrt(stat+syst)<<endl;
	};
      
      // Performs the stability check
      grace_file_t stabilityPlot(baseOut+"/"+ensName+"/stab"+to_string(EStarInGeV)+".xmg");
      for(double statOverSyst=400;statOverSyst>=1e-4;statOverSyst/=2)
	{
	  reconstruct(statOverSyst);
	  stabilityPlot.write_ave_err(log(statOverSyst),R.ave_err());
	}
      
      // Final reconstruction
      reconstruct(exp(-1));
      
      /// Print coefficients
      grace_file_t gPlot(baseOut+"/"+ensName+"/g"+to_string(EStarInGeV)+".xmg");
      for(size_t iT=0;iT<nT;iT++)
	gPlot.write_xy(iT,g[iT]);
      
      /// Estimates the cancellation
      // djack_t maxCoeff;
      // maxCoeff=0;
      // {
      // 	djack_t r;
      // 	r=0;
      // 	for(size_t iT=0;iT<nT;iT++)
      // 	  {
      // 	    const djack_t contr=g[iT]*corr[iT]/preco[iT];
      // 	    r+=contr;
      // 	    maxCoeff=max(abs(contr),maxCoeff);
      // 	  }
      // 	const djack_t cancel=maxCoeff/r;
      // 	cout<<EStar<<" "<<cancel.ave_err()<<endl;
      // }
      // double errEst=0;
      // for(size_t iT=0;iT<nT;iT++)
      // 	for(size_t iS=0;iS<nT;iS++)
      // 	  errEst+=gRes[iT].ave()*gRes[iS].ave()*v[iS+nT*iT]*(iS==iT || 1)/(preco[iT]*preco[iS]);
      // cout.precision(16);
      // cout<<R.err()<<" "<<sqrt(errEst)<<endl;
      
      // Plots the reconstruction of the target function
      grace_file_t recoPlot(baseOut+"/"+ensName+"/reco"+to_string(EStarInGeV)+".xmg");
      recoPlot.write_line([&](const double& EInGeV)
      {
	const double E=EInGeV/aInv.ave();
	double s;
	s=0;
	// double kMin=1e300,kMax=0;
	for(size_t iT=0;iT<nT;iT++)
	  {
	    const double k=b(iT,E)/preco[iT];
	    s+=g[iT]*k;
	    
	    // kMin=min(kMin,k);
	    // kMax=max(kMax,k);
	  }
	s*=sqr(E);
	
	// cout<<E<<" "<<kMin<<" "<<kMax<<endl;
	
	return s;
      },EMin*aInv.ave(),4);
      
      recoPlot.write_line([&](const double& EInGeV)
      {
	const double E=EInGeV/aInv.ave();
	return targ(E)*sqr(E);
      },EMin*aInv.ave(),4);
      
      RofE.emplace_back(EStarInGeV,R);
      
      EStar+=sigma/10;
    }
  
  glbRplot.set_legend(ensName);
  
  for(grace_file_t* plot : {&RPlot,&glbRplot})
    plot->write_polygon([&RofE](const double& E)
  {
    vector<pair<double,size_t>> distPos(RofE.size());
    for(size_t iE=0;iE<RofE.size();iE++)
      distPos[iE]={fabs(get<0>(RofE[iE])-E),iE};
    sort(distPos.begin(),distPos.end());
    
    vector<double> x(3);
    djvec_t y(3);
    for(size_t iE=0;iE<3;iE++)
      tie(x[iE],y[iE])=RofE[distPos[iE].second];
    
    const djvec_t pars=poly_fit(x,y,2);
    
    return poly_eval(pars,E);
  },get<0>(RofE.front()),get<0>(RofE.back()));
  
  RPlot.write_line([c=norm[tNorm].ave()](const double& x)
  {
    return c;
  },EMin,4,grace::BLACK);
}

int main()
{
  set_njacks(50);
  
  grace_file_t RPlot(baseOut+"/R.xmg");
  
  for(const char* ensName : {"B.72.64","B.72.96","C.06.80","C.06.112","D.54.96","E.44.112"})
    analyzeEns(RPlot,ensName);
  
  return 0;
}
