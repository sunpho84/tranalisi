#include <tranalisi.hpp>

#include <gsl/gsl_sf_bessel.h>

constexpr double mpiPhys=0.1350,fpiPhys=0.1304,xiPhys=sqr(mpiPhys/(4.0*M_PI*fpiPhys)),XPhys=0.13406719810514604613;
constexpr bool logFit=true;

double G1t(double x)
{
  const int nterm=20;
  const int mul[nterm]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
  double G1t=0;
  
  for(int nn=0;nn<nterm;nn++)
    {
      double y=x*sqrt(nn+1);
      G1t+=4*mul[nn]*gsl_sf_bessel_K1(y)/y;
      // const double appr=exp(-y)/y*sqrt(M_PI/2);
      // cout<<y<<" "<<gsl_sf_bessel_K1(y)<<" "<<appr<<" "<<gsl_sf_bessel_K1(y)/appr<<endl;
    }
  
  return G1t;
}

  const double dl1b=-0.4,dl2b=4.3,dl3b=2.9,dl4b=4.4;
  
  void volume(double xl,double & g1t,double & g2tm,double & g2tf,double & g2tv)
  {
    double coef1=-55.0/18+4*dl1b+8*dl2b/3-2.5*dl3b-2*dl4b;
    double coef2=112.0/9-8*dl1b/3-32*dl2b/3;
    double coef3=-7.0/9+2*dl1b+4*dl2b/3-3*dl4b;
    double g0=2-0.5*M_PI;
    double g1=0.25*M_PI-0.5;
    double g2=0.5-0.125*M_PI;
    double g3=3.0*M_PI/16-0.5;
    double dmn,b0,b2;
    
    const int nterm=20;
    const int mul[nterm]={6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24};
    int nn;
    double y;
    
    g1t=g2tm=g2tf=g2tv=0;
    
    for(nn=0;nn<nterm;nn++)
      {
	y=xl*sqrt((double)(nn+1));
	dmn=mul[nn];
	b2=gsl_sf_bessel_K0(y);
	g2tv+=dmn*b2;
	b0=gsl_sf_bessel_K1(y);
	b2=(b2+2*b0/y)/y;
	g1t+=4*dmn*b0/y;
	g2tm+=dmn*(coef1*b0+coef2*b2+13*g0*b0/3-(40*g0+32*g1+26*g2)*b2/3)/y;
	g2tf+=2*dmn*(coef3*b0+coef2*b2+(8*g0-13*g1)*b0/6-(40*g0-12*g1-8*g2-13*g3)*b2/3)/y;
      }
  }
  
void corr_FSE_CDH(double &corrvm,double &corrvf,double ampi,double afpi,double lato)
{
  double dn=(4*M_PI*4*M_PI);
      
  double corr0;
  double ampi2=ampi*ampi;
  double xl=lato*ampi;
  double csi=2*ampi2/(dn*afpi*afpi);
  
  double g1t,g2tm,g2tf,g2tv;
  
  volume(xl,g1t,g2tm,g2tf,g2tv);
  
  corr0=-2*(dl4b*g1t-0.5*dl3b*g2tv);
  
  corrvm=1+0.25*csi*g1t-csi*csi*(g2tm-0.25*corr0);
  corrvf=1-csi*g1t+csi*csi*(g2tf-corr0);
}

struct perens_t
{
  string path;
  size_t T;
  size_t L;
  size_t iB;
  size_t TH;
  double am;
  size_t tMin;
  size_t tMax;
  
  djack_t aM,Z2,af,X,xi;
  
  void process()
  {
    cout<<path<<" ";
    const djvec_t corrP5P5=read_djvec(path+"/jacks/mes_contr_P5P5",T).symmetrized();
    
    two_pts_fit(Z2,aM,corrP5P5,TH,tMin,tMax,path+"/plots/effMassP5P5.xmg");
    
    af=sqrt(Z2)*2*am/(aM*sinh(aM));
    
    // cout<<aM.ave_err()<<endl;
    // cout<<af.ave_err()<<endl;
    
    xi=sqr(aM/(4.0*M_PI*af));
    X=pow(pow(aM,4.0)*af,0.2);
    cout<<smart_print(xi.ave_err())<<" "<<smart_print(af.ave_err())<<" "<<smart_print(afCorr().ave_err())<<endl;
  }
  
  
  auto fse_MF() const
  {
    const djack_t xlam=aM*L;
    const djack_t xil=aM*aM/sqr(4*M_PI*af);
    
    djack_t g1t=jackCall(G1t,xlam);
    
    const djack_t corr_mpi=sqrt(1+xil*g1t);
    const djack_t corr_fpi=1-2*xil*g1t;
    
    djack_t corrvm,corrvf;
    jackCall(corr_FSE_CDH,corrvm,corrvf,aM,af,L);
    
    // cout<<corr_mpi<<" "<<corrvm<<endl;
    // cout<<corr_fpi<<" "<<corrvf<<endl;
    
    //return std::make_pair(corr_mpi,corr_fpi);
    return std::make_pair(corrvm,corrvf);
  }
  
  djack_t xiCorr() const
  {
    const auto [fseM,fseF]=fse_MF();
    
    return sqr(aM/fseM/(4.0*M_PI*af/fseF));
  }
  
  djack_t afCorr() const
  {
    const auto [fseM,fseF]=fse_MF();
    
    return af/fseF;
  }
  
  djack_t XCorr() const
  {
    const auto [fseM,fseF]=fse_MF();
    
    return pow(pow(aM/fseM,4.0)*af/fseF,0.2);
  }
  
  perens_t(const string& path) :
    path(path)
  {
    input_file_t input(path+"/input");
    
    T=input.read<size_t>("T");
    TH=L=T/2;
    iB=input.read<size_t>("iB");
    am=input.read<double>("am");
    tMin=input.read<size_t>("tInt");
    tMax=input.read<size_t>();
  }
};

int main()
{
  set_njacks(30);
  
  cout<<"Xi phys: "<<xiPhys<<endl;
  
  vector<perens_t> ensList;
  for(const char* str : {"cA.30.32","cA.40.24","cA.53.24",
			 "cB.14.64","cB.25.24","cB.25.32","cB.25.48","cB.72.64","cB.72.96",
			 "cC.06.80",
			 "cD.54.96"})
    {
      ensList.push_back((string)str);
      ensList.back().process();
    }
  
  auto getEns=[&ensList](const string& name) ->perens_t&
  {
    perens_t* ref=nullptr;
    for(auto& ens : ensList)
      if(ens.path==name)
	{
	  if(ref==nullptr)
	    ref=&ens;
	  else
	    CRASH("Argh!");
	}
    
    return *ref;
  };
  
  const size_t nB=4;
  
  grace_file_t plotXFunXi("plots/XFunXi.xmg");
  grace_file_t plotFPiFunXi("plots/fPiFunXi.xmg");
  grace_file_t plotMPiLFunXi("plots/MPiLFunnXi.xmg");
  grace_file_t plotFPiFunXiCorr("plots/fPiFunXiCorr.xmg");
  grace_file_t plotfPiFunMPiLB1("plots/fPiFunMPiLB1.xmg");

  std::vector<grace_file_t*> plotList{&plotXFunXi,&plotFPiFunXi,&plotMPiLFunXi,&plotFPiFunXiCorr};
  for(auto& plot : plotList)
    plot->set_color_scheme({grace::RED,grace::BLACK,grace::BLUE,grace::GREEN4});
  
  for(size_t iB=0;iB<nB;iB++)
    {
      for(const auto& ens : ensList)
	if(ens.iB==iB)
	  {
	    plotXFunXi.write_ave_err((ens.xiCorr()/xiPhys).ave(),ens.X.ave_err());
	    plotFPiFunXi.write_ave_err((ens.xi/xiPhys).ave(),ens.af.ave_err());
	    plotMPiLFunXi.write_ave_err((ens.xi/xiPhys).ave(),(ens.aM*ens.L).ave_err());
	  }
      
      for(auto& plot : plotList)
	{
	  plot->set_no_line();
	  plot->new_data_set();
	}
    }
  
  for(size_t withWithoutFSE=0;withWithoutFSE<2;withWithoutFSE++)
    {
      const auto y=[withWithoutFSE](const perens_t& ens)
      {
	const auto [fse_M,fse_f]=ens.fse_MF();
	
	if constexpr(0)
	  {
	    if(withWithoutFSE==0)
	      return ens.af;
	    else
	      return ens.af/fse_f;
	  }
	else {
	  if(withWithoutFSE==0)
	    return ens.aM;
	  else
	    return ens.aM/fse_M;
	}
      };
      
      plotfPiFunMPiLB1.set_all_colors(withWithoutFSE?(grace::RED):(grace::BLACK));
      for(const auto& ens : ensList)
	if(ens.iB==1 and ens.am==0.0025)
	  plotfPiFunMPiLB1.write_ave_err(exp(-ens.aM.ave()*ens.L),y(ens).ave_err());
      plotfPiFunMPiLB1.new_data_set();
      
      plotfPiFunMPiLB1.set_all_colors(withWithoutFSE?(grace::RED):(grace::BLACK));
      for(const auto& ens : ensList)
	if(ens.iB==1 and ens.am==0.00072)
	  plotfPiFunMPiLB1.write_ave_err(exp(-ens.aM.ave()*ens.L),y(ens).ave_err());
      plotfPiFunMPiLB1.new_data_set();
    }
  
  for(size_t iB=0;iB<nB;iB++)
    {
      for(const auto& ens : ensList)
	if(ens.iB==iB)
	  plotFPiFunXiCorr.write_ave_err(ens.xiCorr().ave()/xiPhys,ens.afCorr().ave_err());
      
      plotFPiFunXiCorr.set_no_line();
      plotFPiFunXiCorr.new_data_set();
    }
  
  /////////////////////////////////////////////////////////////////
  
  const vector<string> pathFitList={"cB.72.96","cB.14.64","cB.25.48"};
  const vector<string> pathFitList2={"cA.30.32","cA.40.24","cA.53.24"};
  
  auto ansatz=[](const auto& pars,const auto& x)
  {
    if constexpr(logFit)
      {
	const auto xi=(x+1)*xiPhys;
	return pars[0]*(1-2*xi*log(xi/xiPhys)+pars[1]*(xi-xiPhys));
      }
    else
      return pars[0]*(1+pars[1]*x+pars[2]*x*x);
  };
  
  auto fitFSE=[&getEns,&ansatz](const string& plotPath,
				const vector<string>& pathFitList,
				const double min,const double max)
  {
    jack_fit_t fit;
    djvec_t pars(logFit?2:3);
    fit.add_fit_par(pars[0],"A",{0.15,0.01});
    fit.add_fit_par(pars[1],"B",{0.1,0.01});
    if constexpr(not logFit)
      fit.add_fit_par(pars[2],"C",{0.2,0.01});
    
    for(const auto& path : pathFitList)
      {
	const perens_t& ens=getEns(path);
	fit.add_point(ens.afCorr(),[&ansatz,x=ens.xiCorr().ave()/xiPhys-1](const vector<double>& pars,const size_t iel)
	{
	  return ansatz(pars,x);
	});
      }
    
    fit.fit();
    
    grace_file_t fitPlot(plotPath);
    for(const auto& path : pathFitList)
      {
	const perens_t& ens=getEns(path);
	fitPlot.write_ave_err(ens.xiCorr().ave()/xiPhys,ens.afCorr().ave_err());
      }
    
    fitPlot.write_polygon([&pars,&ansatz](double x)
    {
      x-=1;
      return ansatz(pars,x);
    },min,max);
    
    return pars;
  };
  
  const djvec_t pars=fitFSE("plots/fitFPi_vs_Xi.xmg",pathFitList,1,3);
  const djvec_t pars2=fitFSE("plots/fit2FPi_vs_Xi.xmg",pathFitList2,1,5);
  
  cout<<pars.ave_err()<<endl;
  cout<<pars2.ave_err()<<endl;
  
  //for(const auto& path : {"cB.72.64","cB.72.96","cB.14.64","cB.25.32","cB.25.48","cC.06.80","cD.54.96"})
  const vector<ave_err_t> knownA{{0,0},{0.4132,0.0015},{0.3515,0.0013},{0.2921,0.0010}};
  for(const auto& path : {"cB.72.64","cB.72.96","cC.06.80","cD.54.96"})
    {
      const perens_t& ens=getEns(path);
      
      const double x=ens.xiCorr().ave()/xiPhys-1;
      const djack_t retuningFact=ansatz(pars,x)/ansatz(pars,0);
      const djack_t a=ens.afCorr()/fpiPhys/retuningFact;
      
      cout<<"a("<<path<<"): "<<(a.ave_err())<<
	" , mass retuning factor: "<<smart_print(retuningFact.ave_err())<<
	" , FVE correction factor: "<<smart_print((ens.afCorr()/ens.af).ave_err())<<
	" , uncorrected: "<<smart_print((ens.af/fpiPhys).ave_err())<<
	" , rat with known: "<<knownA[ens.iB].ave()/a.ave()<<" "<<knownA[ens.iB].err()/a.ave()<<endl;
    }
  
  return 0;
}
