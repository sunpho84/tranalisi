#include <tantaloBacco.hpp>

#include <tranalisi.hpp>

int main()
{
  set_njacks(15);
  const double lambda=0.0;
  const int hasBwSignal=0;
  // constexpr bool dontUseTantalo=false;
  // constexpr bool useTantalo=true;
  // double sigma=0.2;
  const size_t T=245;
  const size_t tMin=1;
  const double E0=0.00001;
  djvec_t corr(T/2+1);
  corr=0;
  
  PrecFloat::setDefaultPrecision(1024);
  
  using namespace Bacco;
  
  CorrelatorPars correlatorPars(T,hasBwSignal,E0,corr);
  
  const int tMax=128;
  const double Estar=0.5;
  const double sigma=0.25;
  std::unique_ptr<ReconstructionEngine> reconstructor=make_unique<GaussReconstructor>(correlatorPars,tMin,tMax,Estar,lambda,E0,sigma);
  TargettedReconstructor *targettedReconstructor=static_cast<TargettedReconstructor*>(&*reconstructor);
  const Reconstruction reco=
    reconstructor->getReco();
  
  cout<<"Mean: "<<reco.mean().get()<<endl;
  cout<<"WidthOfSquare: "<<reco.widthOfSquare().get()<<endl;
  cout<<"Width: "<<reco.width().get()<<endl;
  cout<<"Square norm: "<<targettedReconstructor->squareNorm().get()<<endl;
  cout<<"Deviation: "<<targettedReconstructor->deviation(reco)<<endl;
  
  grace_file_t RecPlot("/tmp/RecoDelta"+to_string(Estar)+".xmg");
  RecPlot.set_xaxis_label("E");
  RecPlot.set_yaxis_label("\xD");
  RecPlot.set_no_symbol();
  for(double E=0.001;E<=1;E+=0.005)
    RecPlot.write_xy(E,reco.smearingFunction(E).get());
  
  grace_file_t TargPlot("/tmp/TargetDelta"+to_string(Estar)+".xmg");
  TargPlot.set_xaxis_label("E");
  TargPlot.set_yaxis_label("\xD");
  TargPlot.set_no_symbol();
  for(double E=0.001;E<=3;E+=0.0005)
    TargPlot.write_xy(E,targettedReconstructor->targetFunction(E));
  
  // auto get=[&](const size_t& nT,const double& Estar,const bool tantaloUsage)
  // {
  //   TantaloBaccoPars pars(tantaloUsage,T,tMin,nT,E0,lambda,sigma,useBw);
  //   TantaloBaccoRecoEngine recoEngine(pars,Estar);
  //   TantaloBaccoReco reco(recoEngine,Estar,corr);
    
  //   PrecFloat width2=reco.A.formWith(reco.g,reco.g);
  //     PrecFloat ave=0.0;
  //     for(size_t i=0;i<nT;i++)
  // 	ave+=reco.g[i]*reco.M[i];
  //   const PrecFloat width=sqrt(width2);
  //   // for(size_t it=0;it<nT;it++)
  //   //   cout<<reco.g[it].get()<<endl;
    
  //   grace_file_t D("/tmp/RecoDelta"+to_string(Estar)+"_"+to_string(nT)+(tantaloUsage?"use":"dontuse")+".xmg");
  //   D.set_xaxis_label("E");
  //   D.set_yaxis_label("\xD");
  //   D.set_no_symbol();
  //   D.set_subtitle("Nt="+to_string(nT));
  //   for(double E=E0.get();E<=1;E+=0.005)
  //     D.write_xy(E,reco.recoDelta(E).get());
    
  //   if(tantaloUsage)
  //     {
  // 	grace_file_t D("/tmp/TargetDelta"+to_string(Estar)+"_"+to_string(nT)+".xmg");
  // 	D.set_xaxis_label("E");
  // 	D.set_yaxis_label("\xD");
  // 	D.set_no_symbol();
  // 	D.set_subtitle("Nt="+to_string(nT));
  // 	for(double E=E0.get();E<=1;E+=0.005)
  // 	  D.write_xy(E,reco.Delta(Estar,E).get());
	
  // 	grace_file_t Er("/tmp/ErrDelta"+to_string(Estar)+"_"+to_string(nT)+".xmg");
  // 	Er.set_xaxis_label("E");
  // 	Er.set_yaxis_label("\xD");
  // 	Er.set_no_symbol();
  // 	Er.set_subtitle("Nt="+to_string(nT));
  // 	for(double E=E0.get();E<=1;E+=0.005)
  // 	  Er.write_xy(E,(reco.Delta(Estar,E).get()-reco.recoDelta(E)).get());
  //     }
    // cout<<endl;
    // PrecFloat max=0.0;
    // double mP=0;
    // for(size_t iT=0;iT<nT;iT++)
    //   {
    // 	const PrecFloat tmp=abs(reco.g[iT]);
    // 	if(max<tmp)
    // 	  {
    // 	    max=tmp;
    // 	    mP=(double)iT/nT;
    // 	  }
    //   }
    // return std::make_tuple(max,mP);
  //   return width.get();
  // };
  
  // const double w=get(128,0.5,dontUseTantalo);
  // cout<<"w: "<<w<<endl;
  // sigma=w;
  // cout<<sigma<<endl;
  // const double w2=get(128,0.5,useTantalo);
  // cout<<"w2: "<<w2<<endl;
  
  
  // get(16,0.5);
  // get(32,0.5);
  // get(64,0.5);
  
  // return 0;
  
  
  // int prec=64;
  // grace_file_t widthFileVsNt("widthVsNt.xmg");
  // for(size_t nT=8;nT<=1024;nT*=2)
  //   {
  //     cout<<nT<<" "<<prec<<endl;
      
  //     PrecFloat::setDefaultPrecision(prec);
  //     prec*=2;
  //     widthFileVsNt.write_xy(nT,get(nT,0.3));
  //   }
  
  // return 0;
  
  // for(auto nT : {8,16,32,64,128})
  //   {
  //     auto [max,mP]=get(nT,0.4);
  //     cout<<nT<<" "<<max.get()<<endl;
  //   }
  
  // return 0;
  
  // grace_file_t magnitudeFileVsE("magnitudeVsE.xmg");
  // grace_file_t maxPosFileVsE("maxPosVsE.xmg");
  // auto r=[&](const size_t nT)
  //   {
  //     for(double Estar=5;Estar>=0.01;Estar/=1.2)
  // 	{
  // 	  auto [max,mP]=get(nT,Estar);
  // 	  magnitudeFileVsE.write_xy(Estar,max.get());
  // 	  maxPosFileVsE.write_xy(Estar,mP);
  // 	}
      
  //     magnitudeFileVsE.new_data_set();
  //     maxPosFileVsE.new_data_set();
  //   };
  
  // for(auto nT : {8,16,32,64,128})
  //   r(nT);
  
  // PrecFloat sigma=0.025;
  // grace_file_t Dens("/tmp/Dens.xmg");
  
  
  // vector<jack_t<PrecFloat>> corr(T/2+1);
  // ifstream input("/tmp/i");
  // for(size_t t=0;t<=T/2;t++)
  //   {
  //     //const PrecFloat f=exp(-(PrecFloat)t*0.4)+exp(-((PrecFloat)T-t)*0.4);
  //     double f;
  //     input>>f;
  //     corr[t]=f;
  //   }
  // gen_t g(31241);
  // do
  //   {
  //     
  //     const string Etag=to_string(Estar.get());
  //     grace_file_t RecDelta("/tmp/RecDelta"+Etag+".xmg");
  //     grace_file_t RecDelta2("/tmp/RecDelta2"+Etag+".xmg");
  //     grace_file_t TargDelta("/tmp/TargDelta"+Etag+".xmg");
  //     grace_file_t ErrDelta("/tmp/ErrDelta"+Etag+".xmg");
      
  //     PrecFloat y=0;
  //     for(size_t it=0;it<nT;it++)
  // 	y+=reco.g[it]*exp(-0.4*(PrecFloat)it);
      
  //     for(PrecFloat E=0.0;E<Estar*3;E+=0.01)
  // 	{
  // 	  const PrecFloat sTarg=recoEngine.Delta(Estar,E);
  // 	  const PrecFloat sReco=reco.recoDelta(E);
  // 	  RecDelta.write_xy(E.get(),sReco.get());
  // 	  TargDelta.write_xy(E.get(),sTarg.get());
  // 	  const PrecFloat dS=sReco-sTarg;
  // 	  ErrDelta.write_xy(E.get(),dS.get());
	  
  // 	  PrecFloat sReco2=0.0;
  // 	  for(size_t iT=0;iT<nT;iT++)
  // 	    sReco2+=reco.g(iT)*reco.bT(iT+tMin,E);
  // 	  sReco2+=(reco.g(2)*g.get_double(-1e-5,1e-5))*reco.bT(2+tMin,E);
  // 	  RecDelta2.write_xy(E.get(),sReco2.get());
  // 	}
      
  //     const PrecFloat rd=reco.recoDensity()[0];
  //     for(PrecFloat E=Estar-sigma;E<Estar+sigma;E+=0.01)
  // 	Dens.write_xy(E.get(),rd.get());
      
  //     grace_file_t gFile("/tmp/g"+Etag+".xmg");
  //     for(size_t it=0;it<nT;it++)
  // 	gFile.write_xy(it,reco.g[it].get());
      
  //     const PrecFloat EstarPrime=sqr(Estar+sigma)/Estar;
  //     const PrecFloat sigmaPrime=sigma*sqrt(EstarPrime/Estar);
      
  //     //cout<<(sigma+sigmaPrime).get()<<" "<<(EstarPrime-Estar).get()<<endl;
      
  //     Estar=EstarPrime;
  //     sigma=sigmaPrime;
      
  //     cout<<Estar.get()<<" "<<sigma.get()<<endl;
  //   }
  // while(Estar<1);
  
  return 0;
}
