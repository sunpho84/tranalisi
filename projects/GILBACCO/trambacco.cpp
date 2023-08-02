#include <tranalisi.hpp>

using Real=double;
using Jack=jack_t<Real>;
using Jvec=vmeas_t<Jack>;
using Vect=Vector<Real,Dynamic>;
using Matr=Matrix<Real,Dynamic,Dynamic>;
constexpr Real tol=1e-8;

const Real N=2;
const size_t tMin=1;
size_t T;
const size_t tNorm=5;
const Real sigmaInGeVEMin=0.200;
constexpr bool scaleSigma=true;
constexpr Real MPiPhys=0.1350,fPiPhys=0.1304;
//double precoEn;

const string baseIn="/home/francesco/QCD/LAVORI/GM3";
const string baseOut="/home/francesco/QCD/LAVORI/TRAMBACCO";

auto loadCorr(const string& ensName,
	      const string& chann)
{
  raw_file_t fin(baseIn+"/"+ensName+"/data/"+chann,"r");
  const size_t nConfs=fin.bin_read<size_t>();
  
  [[maybe_unused]] const size_t nHits=fin.bin_read<size_t>();
  T=fin.bin_read<size_t>();
  
  index_t iRaw({{"nConfs",nConfs},{"T/2+1",T/2+1}});
  vector<double> raw(iRaw.max());
  fin.bin_read(raw);
  
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
  
  // constexpr bool preciseCov=true;
  
  djvec_t covMatr((T/2+1)*(T/2+1));
  covMatr=0;
  
  // if constexpr(not preciseCov)
  //   {
  //     for(size_t t=0;t<=T/2;t++)
  // 	for(size_t s=0;s<=T/2;s++)
  // 	  covMatr[t+(T/2+1)*s]=cov(res[t],res[s]);
  //   }
  // else
  //   {
  gen_t gen(232452);
  const size_t nBoots=1000;
  const double clustSize=nConfs/(double)njacks;
  for(size_t iJack=0;iJack<=njacks;iJack++)
    {
      vector<double> ave(T/2+1,0.0);
      vector<double> c((T/2+1)*(T/2+1));
      vector<double> b(T/2+1,0.0);
      
      cout<<iJack<<" [0,"<<iJack*clustSize<<") ["<<(iJack+1)*clustSize<<" "<<nConfs<<")"<<endl;
      cout<<"......"<<endl;
      
      for(auto& ci : c)
	ci=0.0;
      
      for(size_t iBoot=0;iBoot<nBoots;iBoot++)
	{
	  for(auto& bi : b)
	    bi=0.0;
	  
	  const size_t nRedConfs=nConfs-(iJack!=njacks)*clustSize;
	  for(size_t iConf=0;iConf<nRedConfs;iConf++)
	    {
	      const double i=gen.get_double(0,nConfs-clustSize);
	      const size_t jConf=i+(i>clustSize*iJack)*clustSize;
	      for(size_t t=0;t<=T/2;t++)
		b[t]+=raw[iRaw({jConf,t})];
	    }
	  
	  for(size_t t=0;t<=T/2;t++)
	    b[t]/=nRedConfs;
	  
	  for(size_t t=0;t<=T/2;t++)
	    {
	      ave[t]+=b[t];
	      for(size_t s=0;s<=T/2;s++)
		c[s+(T/2+1)*t]+=b[t]*b[s];
	    }
	}
      
      for(size_t t=0;t<=T/2;t++)
	{
	  ave[t]/=nBoots;
	  for(size_t s=0;s<=T/2;s++)
	    c[s+(T/2+1)*t]/=nBoots;
	}
      
      for(size_t t=0;t<=T/2;t++)
	for(size_t s=0;s<=T/2;s++)
	  c[t+(T/2+1)*s]-=ave[s]*ave[t];
      
      for(size_t t=0;t<=T/2;t++)
	{
	  for(size_t s=t+1;s<=T/2;s++)
	    {
	      c[t+(T/2+1)*s]/=sqrt(c[t+(T/2+1)*t])*sqrt(c[s+(T/2+1)*s]);
	      c[s+(T/2+1)*t]=c[t+(T/2+1)*s];
	    }
	  c[t+(T/2+1)*t]=1;
	}
      
      for(size_t t=0;t<=T/2;t++)
	for(size_t s=0;s<=T/2;s++)
	  covMatr[t+(T/2+1)*s][iJack]=c[t+(T/2+1)*s];
      // for(size_t t=0;t<=T/2-10;t++)
      // 	{
      // 	  for(size_t d=0;d<=10;d++)
      // 	    {
      // 	      const size_t s=(t+d)%(T/2+1);
      // 	      cout<<pow(covMatr[t+(T/2+1)*s],1.0/(std::max(1,(int)d)))<<" ";
      // 	    }
      // 	  cout<<endl;
      // 	}
    }
  
  // covMatr.clusterize(clustSize);
  
  const size_t dMax=10;
  djvec_t rs(dMax);
  for(size_t d=0;d<dMax;d++)
    {
      djack_t su;
      su=0;
      size_t n=0;
      for(size_t t=1;t<=T/2-dMax;t++)
	{
	  const size_t s=(t+d)%(T/2+1);
	  su+=covMatr[t+(T/2+1)*s];
	  n++;
	}
      su/=n;
      
      rs[d]=su;
    }
  rs.ave_err().write("/tmp/fitto.xmg");
  cout<<endl;
  
  for(size_t t=0;t<=T/2;t++)
    for(size_t s=0;s<=T/2;s++)
      {
	const size_t d=abs((int)t-(int)s);
	if(d<dMax) covMatr[t+(T/2+1)*s]=rs[d];
	else covMatr[t+(T/2+1)*s]=0.0;
	
    }

  for(size_t t=2;t<=T/2-10;t++)
    {
      for(size_t d=0;d<=10;d++)
	{
	  const size_t s=(t+d)%(T/2+1);
	  cout<<smart_print(covMatr[t+(T/2+1)*s].ave_err())<<" ";
	}
      cout<<endl;
    }
  cout<<endl;
  
  for(size_t t=0;t<=T/2;t++)
    {
      const double d=res[t].err();
      
      for(size_t s=0;s<=T/2;s++)
	{
	  covMatr[t+(T/2+1)*s]*=d;
	  covMatr[s+(T/2+1)*t]*=d;
	}
    }
  
  return make_tuple(res,covMatr);
}

constexpr int fact(const int& n)
{
  int r=1;
  
  for(int i=2;i<=n;i++)
    r*=i;
  
  return r;
}

Real gauss(const Real& x0,
	   const Real& sigma,
	   const Real& x)
{
  return exp(-sqr((x-x0)/sigma)/2)/(sqrt(2*M_PI)*sigma);
}

Real b(const size_t& it,
       const Real& E)
{
  const Real t=it+tMin;
  
  return exp(-t*E)+exp(-(T-t)*E);
}

void analyzeEns(grace_file_t& glbRplot,
		const string& ensName)
{
  // Don't move, needed to set T
  const auto [cP5P5,dum]=
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
  
  const auto [cVKVK,covVKVK]=
    loadCorr(ensName,"ll_TM_VKVK");
  CRASH("");
  effective_mass(cVKVK).ave_err().write(baseOut+"/"+ensName+"/rho.xmg");
  Jvec corr(nT);
  vector<Real> corrCov(nT*nT);
  // precoEn=effective_mass(cVKVK)[15].ave();
  for(size_t iT=0;iT<nT;iT++)
    {
      const Real fact=sqr(z);//*12*sqr(M_PI)*4/9
      corr[iT]=cVKVK[iT+tMin]*fact;////b(iT,precoEn);
      for(size_t iS=0;iS<nT;iS++)
	corrCov[iT+nT*iS]=covVKVK[(iT+tMin)+(T/2+1)*(iS+tMin)].ave()*sqr(fact);///(b(iT,precoEn)*b(iS,precoEn));
    }
  corr.ave_err().write(baseOut+"/"+ensName+"/rhoPreco.xmg");
  
  vector<Real> preco(nT);
  for(size_t iT=0;iT<nT;iT++)
    preco[iT]=corr[iT].ave();
  
  grace_file_t RPlot(baseOut+"/"+ensName+"/R.xmg");
  const Real EMin=2*aMPi.ave(),EMax=4/aInv.ave();
  
  Jvec norm(nT);
  for(size_t iT=0;iT<nT;iT++)
    norm[iT]=corr[iT]/gslIntegrateUpToInfinity([iT](const Real& E)
    {
      return pow(E,N)*b(iT,E);///b(iT,precoEn);
    },EMin,tol);
  
  norm.ave_err().write(baseOut+"/"+ensName+"/norm.xmg");
  
  grace_file_t sigmaPlot(baseOut+"/"+ensName+"/sigma.xmg");
  
  Real EStar=EMin;
  vector<tuple<Real,Jack>> RofE;
  while(EStar<EMax)
    {
      const Real EStarInGeV=EStar*aInv.ave();
      const Real sigma=sigmaInGeVEMin/aInv.ave()*(scaleSigma?sqrt(EStar/EMin):1);
      const Real sigmaInGeV=sigma*aInv.ave();
      cout<<"E: "<<EStarInGeV<<" GeV, Sigma: "<<sigmaInGeV<<" GeV"<<endl;
      
      sigmaPlot.write_xy(EStarInGeV,sigmaInGeV);
      
      /// Target function
      const auto targ=
	[EStar,sigma](const Real& E)
	{
	  return gauss(EStar,sigma,E)/sqr(E);
	};
      
      // const Real D=
      // 	gslIntegrateUpToInfinity([targ](const Real& E)
      // 	{
      // 	  return sqr(targ(E))*pow(E,2*N);
      // 	},EMin,tol);
      
      Vect e(nT);
      Matr f(nT,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  e(iT)=
	    gslIntegrateUpToInfinity([targ,iT](const double& E)
	    {
	      return targ(E)*pow(E,2*N)*b(iT,E);///b(iT,precoEn);
	    },EMin,tol);
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      gslIntegrateUpToInfinity([iT,iS](const double& E)
	      {
		return pow(E,2*N)*b(iT,E)*b(iS,E);///(b(iT,precoEn)*b(iS,precoEn));
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
	      constexpr Real corrFading=0.9;
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
	  
	  // const Real syst=g.transpose()*QSyst*g;
	  // const Real stat=(Real)(g.transpose()*QStat*g)/statOverSyst;
	  
	  // cout<<statOverSyst<<" "<<sqrt(stat+syst)<<endl;
	};
      
      // Performs the stability check
      grace_file_t stabilityPlot(baseOut+"/"+ensName+"/stab"+to_string(EStarInGeV)+".xmg");
      for(Real statOverSyst=400;statOverSyst>=1e-4;statOverSyst/=2)
	{
	  reconstruct(statOverSyst);
	  stabilityPlot.write_ave_err(log(statOverSyst),R.ave_err());
	}
      
      // ofstream(baseOut+"/"+ensName+"/QSyst"+to_string(EStarInGeV)+".xmg")<<QSyst<<endl;
      // ofstream(baseOut+"/"+ensName+"/QStat"+to_string(EStarInGeV)+".xmg")<<QStat<<endl;
      // ofstream(baseOut+"/"+ensName+"/L"+to_string(EStarInGeV)+".xmg")<<L<<endl;
      
      // Final reconstruction
      reconstruct(exp(-1));
	  // MatrixXd Q=
	  //   QSyst+QStat/exp(-1);
	  // cerr<<Q.inverse()<<endl<<endl;
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
      
      // Plots the reconstruction of the target function
      grace_file_t recoPlot(baseOut+"/"+ensName+"/reco"+to_string(EStarInGeV)+".xmg");
      recoPlot.write_line([&](const double& EInGeV)
      {
	const Real E=EInGeV/aInv.ave();
	Real s;
	s=0;
	// Real kMin=1e300,kMax=0;
	for(size_t iT=0;iT<nT;iT++)
	  {
	    const Real k=b(iT,E)/preco[iT];///b(iT,precoEn)
	    s+=g[iT]*k;
	    
	    // kMin=min(kMin,k);
	    // kMax=max(kMax,k);
	  }
	s*=sqr(E);
	
	// cout<<E<<" "<<kMin<<" "<<kMax<<endl;
	
	return s;
      },EMin*aInv.ave(),EMax*aInv.ave());
      
      recoPlot.write_line([&](const Real& EInGeV)
      {
	const Real E=EInGeV/aInv.ave();
	return targ(E)*sqr(E);
      },EMin*aInv.ave(),EMax*aInv.ave());
      
      RofE.emplace_back(EStarInGeV,R);
      
      EStar+=sigma/10;
    }
  
  glbRplot.set_legend(ensName);
  
  for(grace_file_t* plot : {&RPlot,&glbRplot})
    plot->write_polygon([&RofE](const Real& E)
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
  },get<0>(RofE.front()),get<0>(RofE.back()));
  
  RPlot.write_line([c=norm[tNorm].ave()](const Real& x)
  {
    return c;
  },EMin*aInv.ave(),EMax*aInv.ave(),grace::BLACK);
}

int main()
{
  
  set_njacks(50);
  
  grace_file_t RPlot(baseOut+"/R.xmg");
  
  for(const char* ensName : {"B.72.64","B.72.96","C.06.80","C.06.112","D.54.96","E.44.112"})
  // // for(const char* ensName : {"C.06.80"})
  // for(const char* ensName : {"B.72.64"})
    analyzeEns(RPlot,ensName);
  
  return 0;
}
