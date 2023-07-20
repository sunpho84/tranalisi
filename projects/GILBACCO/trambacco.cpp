#include <tranalisi.hpp>

const double N=2;
const size_t tMin=1,tMax=30;
const size_t nT=tMax-tMin;
size_t T;
const size_t tNorm=5;
const double sigmaInGeV=0.450;

djvec_t loadCorr(const string& chann)
{
  raw_file_t fin("/home/francesco/QCD/LAVORI/GM3/C.06.112/data/"+chann,"r");
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

const double b(const size_t& it,
	       const double& E)
{
  const double t=it+1;
  return exp(-t*E)+exp(-(T-t)*E);
}

int main()
{
  set_njacks(50);
  
  const djvec_t cP5P5=
    loadCorr("ll_TM_P5P5");
  const djack_t aMPi=
    constant_fit(effective_mass(cP5P5),30,40,"/tmp/pion.xmg");
  cout<<"aMPi: "<<aMPi.ave_err()<<endl;
  
  const djack_t aInv=
    0.135/aMPi;
  const djack_t aInFm=
    0.197/aInv;
  cout<<"a^-1: "<<smart_print(aInv)<<" GeV"<<endl;
  cout<<"a: "<<smart_print(aInFm)<<" fm"<<endl;
  
  const djvec_t cVKVK=
    loadCorr("ll_TM_VKVK");
  effective_mass(cVKVK).ave_err().write("/tmp/rho.xmg");
  djvec_t corr(nT);
  for(size_t iT=0;iT<nT;iT++)
    corr[iT]=cVKVK[iT+tMin];
  
  // grace_file_t plot("/tmp/p.xmg");
  // plot.write_line([](const double& x){return ;},-10,10);
  
  const double sigma=sigmaInGeV/aInv.ave();
  const double E0=2*aMPi.ave();
  cout<<"Sigma: "<<sigma<<endl;
  
  grace_file_t RPlot("/tmp/R.xmg");
  const double EMin=E0,EMax=4/aInv.ave();
  
  djvec_t norm(nT);
  for(size_t iT=0;iT<nT;iT++)
    norm[iT]=corr[iT]/gslIntegrateUpToInfinity([t=iT+tMin](const double& E)
    {
      return pow(E,N)*exp(-E*t);
    },0);
  
  norm.ave_err().write("/tmp/norm.xmg");
  
  for(double EStar=EMin;EStar<EMax;EStar+=0.01)
    {
      const auto targ=
	[EStar,sigma](const double& E)
	{
	  return gauss(EStar,sigma,E)/pow(E,N);
	};
      
      // const double D=
      // 	gslIntegrateUpToInfinity([targ](const double& E)
      // 	{
      // 	  return sqr(targ(E))*pow(E,2*N);
      // 	},E0);
      
      VectorXd e(nT);
      MatrixXd f(nT,nT);
      
      for(size_t iT=0;iT<nT;iT++)
	{
	  e(iT)=
	    gslIntegrateUpToInfinity([targ,iT](const double& E)
	    {
	      return targ(E)*pow(E,2*N)*b(iT,E);
	    },E0);
	  // cout<<iT<<" "<<e(iT)<<endl;
	  
	  for(size_t iS=iT;iS<nT;iS++)
	    f(iT,iS)=f(iS,iT)=
	      gslIntegrateUpToInfinity([iT,iS](const double& E)
	      {
		return pow(E,2*N)*b(iT,E)*b(iS,E);
	      },E0);
	}
      
      vector<double> v(nT);
      for(size_t iT=0;iT<nT;iT++)
	v[iT]=sqr(corr[iT].err());
      
      djack_t R;
      djvec_t gRes(nT);
      for(size_t iJack=0;iJack<=njacks;iJack++)
	{
	  vector<double> c(nT);
	  for(size_t iT=0;iT<nT;iT++)
	    c[iT]=corr[iT][iJack];
	  const double c2=
	    sqr(norm[tNorm].ave());
	  
	  VectorXd L(nT);
	  MatrixXd Q(nT,nT);
	  
	  for(size_t iT=0;iT<nT;iT++)
	    {
	      L(iT)=c2*e(iT)/c[iT];
	      for(size_t iS=0;iS<nT;iS++)
		Q(iT,iS)=c2*f(iT,iS)/(c[iT]*c[iS]);
	      
	      Q(iT,iT)+=v[iT]/sqr(c[iT]);
	    }
	  
	  const VectorXd g=
	    Q.inverse()*L;
	  
	  R[iJack]=0;
	  for(size_t iT=0;iT<nT;iT++)
	    {
	      gRes[iT][iJack]=g[iT];
	      R[iJack]+=g[iT];
	    }
	  
	  //cout<<D<<endl;
	}
      
      RPlot.write_ave_err(EStar*aInv.ave(),R.ave_err());
      
      grace_file_t recoPlot("/tmp/reco"+to_string(EStar)+".xmg");
      recoPlot.write_line([&corr,&g=gRes](const double& E)
      {
	double s;
	s=0;
	for(size_t iT=0;iT<nT;iT++)
	  s+=b(iT,E)*g[iT].ave()/corr[iT].ave();
	s*=pow(E,N);
	
	return s;
      },E0,4);
      
      recoPlot.write_line([targ](const double& E)
      {
	return targ(E)*pow(E,N);
      },E0,4);
    }
  
  RPlot.write_line([c=norm[tNorm].ave()](const double& x)
  {
    return c;
  },E0,4,grace::BLACK);
  
  return 0;
}
