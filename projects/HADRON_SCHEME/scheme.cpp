#include <tranalisi.hpp>

#include <set>

string W0FrADir,afPionDir;

void new_section(const string section)
{
  cout<<"///////////////////////////////// "<<section<<" ////////////////////////////////"<<endl;
}

const size_t nb=3;
const string beta_tag[nb]={"1.90","1.95","2.10"};

/// List of ml ordered by ib
map<size_t,set<double>> amlList;

const bool usePhi=true;

const double MPionZeroExp=0.13498;
const double fPionExp=0.13041;
const double MKaonPlusExp=0.49368;
const double MKaonZeroExp=0.49761;
const double MKaonBarExp=sqrt((sqr(MKaonPlusExp)+sqr(MKaonZeroExp))/2);
const double MEtaExp=sqrt(2*sqr(MKaonBarExp)-sqr(MPionZeroExp));
const double MPhiExp=1.0194;
const double MOmegaExp=1.672;
const double MOmegaFrMPhiExp=MOmegaExp/MPhiExp;
const double MProtonExp=0.938;
const double W0Exp=0.171;
const double msFrMlExp=27.23;

const double MOnlySExp=(usePhi?MPhiExp:MOmegaExp);
const double ratioSExp=sqr(MEtaExp/MOnlySExp);;
const double ratioLExp=sqr(MPionZeroExp/MOnlySExp);
const string ssName=(usePhi?"phi":"omga");

/// Solve for the x of a parabola passing through d
djack_t parab_solve(const djvec_t& pars,const djack_t d,const bool first=true)
{
  const djack_t& a=pars[2],&b=pars[1],&c=pars[0]-d;
  const djack_t& x0=(-b-sqrt(b*b-4*a*c))/(2*a);
  const djack_t& x1=(-b+sqrt(b*b-4*a*c))/(2*a);
  
  if(first)
    return x0;
  else
    return x1;
}

djack_t effMassRemoveExc(const djvec_t& C,const size_t tmin0,const size_t tmax0,const size_t tmin1,const size_t TH,const string& plot_path)
{
  djack_t M0_ext,Z2_0_ext;
  two_pts_fit(Z2_0_ext,M0_ext,C,TH,tmin0,tmax0,plot_path+"_first.xmg");
  djvec_t C1=C;
  for(size_t t=0;t<=TH;t++)
    C1[t]-=two_pts_corr_fun(Z2_0_ext,M0_ext,TH,t,1);
  djack_t M1_ext,Z2_1_ext;
  two_pts_fit(Z2_1_ext,M1_ext,C1,TH,5,tmin0,plot_path+"_second.xmg");
  const djack_t dM10_ext=M1_ext-M0_ext;
  
  jack_fit_t fit;
  djvec_t fit_pars(4);
  fit.add_fit_par_limits(fit_pars[0],"Z2_0",Z2_0_ext.ave(),Z2_0_ext.err(),Z2_0_ext.ave()-4.0*Z2_0_ext.err(),Z2_0_ext.ave()+4*Z2_0_ext.err());
  fit.add_fit_par_limits(fit_pars[1],"M0",M0_ext.ave(),M0_ext.err(),M0_ext.ave()-4*M0_ext.err(),M0_ext.ave()+4*M0_ext.err());
  fit.add_fit_par_limits(fit_pars[2],"Z2_1",Z2_1_ext.ave(),Z2_1_ext.err(),0.0,10*Z2_1_ext.ave());
  fit.add_fit_par_limits(fit_pars[3],"dM10",dM10_ext.ave(),dM10_ext.err(),0.0,10*M0_ext.ave());
  
  for(size_t t=tmin1;t<tmax0;t++)
    fit.add_point(C[t],[t,TH](const vector<double>& p,const size_t ijack)->double
		       {
			 return
			   two_pts_corr_fun(p[0],p[1],TH,t,1)+
			   two_pts_corr_fun(p[2],p[3]+p[1],TH,t,1);
		       });
  
  fit.fit();
  cout<<"combo_fit_pars: "<<endl<<fit_pars.ave_err()<<endl;
  
  /// Plot determine phys for each beta
  grace_file_t plot(plot_path);
  plot.write_vec_ave_err(effective_mass(C).ave_err());
  plot.write_constant_band(tmin0,tmax0,M0_ext);
  plot.write_constant_band(tmin1,tmax0,fit_pars[1]);
  
  return fit_pars[1];
}

/// Esnsemble quantities
struct perens_t
{
  size_t ib;
  string name;
  size_t nq;
  vector<double> amq;
  size_t L,T,TH;
  size_t tminPion,tmaxPion;
  size_t tminKaon,tmaxKaon;
  size_t tminPhi,tmaxPhi;
  size_t tminNucleon,tmaxNucleon;
  size_t tminOmega,tmaxOmega;
  
  perens_t(const string& name) : name(name)
  {
    input_file_t input(name+"/input");
    
    ib=input.read<size_t>("IB");
    
    L=input.read<size_t>("L");
    T=input.read<size_t>("T");
    TH=T/2;
    nq=input.read<int>("amq");
    amq.resize(nq);
    for(size_t iq=0;iq<nq;iq++)
      amq[iq]=input.read<double>();
    
    tminPion=input.read<size_t>("tintPion");
    tmaxPion=input.read<size_t>();
    
    tminKaon=input.read<size_t>("tintKaon");
    tmaxKaon=input.read<size_t>();
    
    tminPhi=input.read<size_t>("tintPhi");
    tmaxPhi=input.read<size_t>();
    
    tminNucleon=input.read<size_t>("tintNucleon");
    tmaxNucleon=input.read<size_t>();
    
    tminOmega=input.read<size_t>("tintOmega");
    tmaxOmega=input.read<size_t>();
  }
  
  /// Load a barionic contraction
  djvec_t load_bar(const string& tag,const int iproj,const int iWick) const
  {
    string path=combine("%s/jacks/bar_alt_contr_%s_proj_%d_Wick_%d",name.c_str(),tag.c_str(),iproj,iWick);
    
    return -read_djvec(path,T).symmetrized(-1);
  }
  
  /// Return the tag to load a contraction
  static string bar_tag(const size_t iq)
  {
    return combine("SM_q%zu_q%zu_q%zu",iq,iq,iq);
  }
  
  /// Loads the nucleon
  djvec_t load_nucleon(const size_t iq) const
  {
    const string tag=bar_tag(iq);
    
    const djvec_t dir_ii=load_bar(tag,1,0);
    const djvec_t exc_ii=load_bar(tag,1,1);
    const djvec_t dir_ij=load_bar(tag,2,0);
    const djvec_t exc_ij=load_bar(tag,2,1);
    
    const djvec_t ii_nuc=dir_ii-exc_ii;
    const djvec_t ij_nuc=dir_ij-exc_ij;
    const djvec_t nucleon=ii_nuc-ij_nuc;
    
    return nucleon;
  }
  
  /// Reads the omega
  djvec_t load_omega(const size_t iq) const
  {
    const string tag=bar_tag(iq);
    
    const djvec_t dir_ii=load_bar(tag,1,0);
    const djvec_t exc_ii=load_bar(tag,1,1);
    const djvec_t dir_ij=load_bar(tag,2,0);
    const djvec_t exc_ij=load_bar(tag,2,1);
    const djvec_t ii=dir_ii-2*exc_ii;
    const djvec_t ij=dir_ij-2*exc_ij;
    const djvec_t omega=ii+0.5*ij;
    
    return omega;
  }
  
  string meson_file_name(const size_t iq1,const size_t iq2,const string& channel) const
  {
    return combine("%s/jacks/mes_contr_%s_SM_q%zu__SM_q%zu",name.c_str(),channel.c_str(),iq1,iq2);
  }
  
  /// Reads a meson
  djvec_t load_meson(const size_t iq1,const size_t iq2,const string& channel) const
  {
    return read_djvec(meson_file_name(iq1,iq2,channel),T).symmetrized();
  }
  
  djack_t w0_fr_a;
  
  djvec_t cPion,ePion;
  djack_t amPion,afPion,aPhiPion;
  
  djvec_t cNucleon,eNucleon;
  djack_t amNucleon;
  
  vector<djvec_t> cKaon,eKaon;
  djvec_t amKaon;
  
  vector<djvec_t> cOmega,eOmega;
  djvec_t amOmega;
  
  vector<djvec_t> cPhi,ePhi;
  djvec_t amPhi;
  
  djvec_t amOmegaFrAmPhi;
  
  djvec_t amOnlyS;
  
  djvec_t amEtaSS;
  
  djvec_t ratioS;
  
  djvec_t ratioL;
  
  djack_t amsPhys;
  
  // djack_t amsPhiPhys;
  
  /// Load the correlators
  perens_t& loadCorrs()
  {
    /// Name of the dir where to put the plots
    const string dir=name+"/plots/corrs";
    
    cPion=load_meson(0,0,"P5P5");
    cPion.ave_err().write(dir+"/pion.xmg");
    cNucleon=load_nucleon(0);
    
    for(auto& c : {&cKaon,&cOmega,&cPhi})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	cOmega[iq]=load_omega(iq);
	cOmega[iq].ave_err().write(combine("%s/omega_%zu.xmg",dir.c_str(),iq));
	cKaon[iq]=load_meson(0,iq,"P5P5");
	cKaon[iq].ave_err().write(combine("%s/kaon_%zu.xmg",dir.c_str(),iq));
      }
    
    for(size_t iq=0;iq<nq;iq++)
      {
	cPhi[iq]=0.0;
	for(size_t i=1;i<=3;i++)
	  cPhi[iq]=load_meson(iq,iq,"V"+to_string(i)+"V"+to_string(i));
	cPhi[iq]/=3.0;
	cPhi[iq].ave_err().write(combine("%s/phi_%zu.xmg",dir.c_str(),iq));
      }
    
    return *this;
  }
  
  /// Load w0/a
  perens_t& loadW0fra()
  {
    static int progressive=12423;
    
    double ave,err;
    ifstream w0_fr_a_file(W0FrADir+"/"+name+"/w0_fr_a.txt");
    w0_fr_a_file>>ave>>err;
    w0_fr_a.fill_gauss(ave,err,progressive++);
    
    cout<<"W0/a: "<<w0_fr_a.ave_err()<<endl;
    
    return *this;
  }
  
  /// Load afPi
  perens_t& loadAfPion()
  {
    static int progressive=7457343;
    
    double ave,err;
    ifstream afPion_file(afPionDir+"/"+name+"/f.txt");
    afPion_file>>ave>>err;
    afPion.fill_gauss(ave,err,progressive++);
    
    cout<<"af: "<<afPion.ave_err()<<endl;
    
    return *this;
  }
  
  /// Computes the effective masses
  perens_t& effectiveMasses()
  {
    ePion=effective_mass(cPion);
    eNucleon=effective_mass(cNucleon,TH,0);
    
    for(auto& c : {&eKaon,&eOmega,&ePhi})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	eOmega[iq]=effective_mass(cOmega[iq],TH,0);
	eKaon[iq]=effective_mass(cKaon[iq]);
	ePhi[iq]=effective_mass(cPhi[iq]);
      }
    
    return *this;
  }
  
  /// Fit the effective masses
  perens_t& fitMasses()
  {
    /// Dir where to put plots
    const string dir=name+"/plots/effMass";
    
    amPion=constant_fit(ePion,tminPion,tmaxPion,dir+"/pion.xmg");
    aPhiPion=afPion*sqrt(amPion);
    amNucleon=constant_fit(eNucleon,tminNucleon,tmaxNucleon,dir+"/nucleon.xmg");
    
    for(auto& c : {&amKaon,&amOmega,&amPhi,&amOmegaFrAmPhi,&amOnlyS})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	const double aApprox[3]={0.0886,0.0815,0.0619};
	const size_t amOmegaFrAmPhiTMin=7*aApprox[0]/aApprox[ib];
	const size_t amOmegaFrAmPhiTMax=11*aApprox[0]/aApprox[ib];
	
	amOmega[iq]=constant_fit(eOmega[iq],tminOmega,tmaxOmega,combine("%s/omega_%zu.xmg",dir.c_str(),iq));
	amKaon[iq]=constant_fit(eKaon[iq],tminKaon,tmaxKaon,combine("%s/kaon_%zu.xmg",dir.c_str(),iq));
	amPhi[iq]=constant_fit(ePhi[iq],tminPhi,tmaxPhi,combine("%s/phi_%zu.xmg",dir.c_str(),iq));
	amOmegaFrAmPhi[iq]=constant_fit((djvec_t)(eOmega[iq]/ePhi[iq]),amOmegaFrAmPhiTMin,amOmegaFrAmPhiTMax,combine("%s/omegaFrPhi_%zu.xmg",dir.c_str(),iq));
	
	if(usePhi)
	  amOnlyS[iq]=amPhi[iq];
	else
	  amOnlyS[iq]=amOmega[iq];
	
	const djack_t openPhi=amPhi[iq]-2*sqrt(sqr(amKaon[iq])+sqr(2*M_PI/L));
	cout<<"Phi decay channel mass difference: "<<openPhi.ave_err()<<endl;
      }
    
    return *this;
  }
  
  /// Print fitted masses
  perens_t& printMasses()
  {
    cout<<"amPion: "<<amPion<<endl;
    cout<<"amNucleon: "<<amNucleon<<endl;
    
    auto print=[this](const string& name,const djvec_t& val)
	       {
		 cout<<name<<":"<<endl;
		 for(size_t iq=1;iq<nq;iq++)
		   cout<<" "<<iq<<": "<<val[iq]<<endl;
	       };
    print("amOmega",amOmega);
    print("amKaon",amKaon);
    
    return *this;
  }
  
  /// Set the ratios
  perens_t& setRatios()
  {
    for(auto& c : {&amEtaSS,&ratioS,&ratioL})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	amEtaSS[iq]=sqrt(2*sqr(amKaon[iq])-sqr(amPion));
    	ratioS[iq]=sqr((djack_t)(amEtaSS[iq]/amOnlyS[iq]));
	ratioL[iq]=sqr((djack_t(amPion/amOnlyS[iq])));
      }
    
    return *this;
  }
  /// Computes ams phys
  perens_t& getAmsPhys()
  {
    // const djvec_t ratioSphiPars=poly_fit(amq,ratioSphi,2,amq[1]*0.9);
    // ratioSphiPlot.write_polygon([ratioSphiPars](const double x){return poly_eval(ratioSphiPars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::GREEN4);
    // const djack_t amsQuadPhi=parab_solve(ratioSphiPars,ratioSphiExp,false);
    // cout<<"ams(quad,phi) "<<name<<": "<<amsQuadPhi.ave_err()<<endl;
    
    // ratioSphiPlot.write_ave_err(amsQuadPhi.ave_err(),{ratioSphiExp,0.0});
    
    // amsPhiPhys=amsQuadPhi;
	
    grace_file_t ratioSPlot(name+"/plots/ratios/s.xmg");
    
    for(size_t iq=0;iq<nq;iq++)
      ratioSPlot.write_ave_err(amq[iq],ratioS[iq].ave_err());
    
    /// Fit the ratioS with a second order polynomial
    const djvec_t ratioSPars=poly_fit(amq,ratioS,2,amq[1]*0.9);
    ratioSPlot.write_polygon([ratioSPars](const double x){return poly_eval(ratioSPars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::GREEN4);
    const djack_t amsQuad=parab_solve(ratioSPars,ratioSExp,false);
    cout<<"ams(quad) "<<name<<": "<<amsQuad.ave_err()<<endl;
    
    // /// Fit the ratioS at all ms with a first order polynomial
    // const djvec_t ratioSLin3Pars=poly_fit(amq,ratioS,1,amq[1]*0.9);
    // ratioSPlot.write_polygon([ratioSLin3Pars](const double x){return poly_eval(ratioSLin3Pars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    // const djack_t amsLin3=-(ratioSLin3Pars[0]-ratioSExp)/ratioSLin3Pars[1];
    // cout<<"ams(lin3): "<<amsLin3.ave_err()<<endl;
    
    // /// Fit the ratioS at all ms with a first order polynomial
    // const djvec_t ratioSLin2Pars=poly_fit(amq,ratioS,1,amq[2]*0.9);
    // ratioSPlot.write_polygon([ratioSLin2Pars](const double x){return poly_eval(ratioSLin2Pars,x);},amq[2]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    // const djack_t amsLin2=-(ratioSLin2Pars[0]-ratioSExp)/ratioSLin2Pars[1];
    // cout<<"ams(lin2): "<<amsLin2.ave_err()<<endl;
    
    amsPhys=amsQuad;
    
    ratioSPlot.write_ave_err(amsPhys.ave_err(),{ratioSExp,0.0});
    
    return *this;
  }
};

vector<perens_t> ensList,ensListUncorr;

double fitEvalX(const double& t)
{
  return t;
}

const vector<double>& fitEvalX(const vector<double>& t)
{
  return t;
}

ave_err_t fitEvalX(const djack_t& t)
{
  return t.ave_err();
}

/// Polynomial fit and evaluation
template <typename T>
djack_t fitAndEval(const string& path,const vector<double>& x,const djvec_t& y,const size_t order,const double min,const double max,const T x0,const string& xlabel="",const string& ylabel="")
{
  const djvec_t pars=poly_fit(x,y,order,min,max);
  const djack_t out=poly_eval(pars,x0);
  
  grace_file_t plot(path);
  plot.write_polygon([&pars](const double& x){return poly_eval(pars,x);},min,max);
  plot.set_no_line();
  plot.write_vec_ave_err(x,y.ave_err());
  plot.set_all_colors(grace::GREEN4);
  plot.write_ave_err(fitEvalX(x0),out.ave_err());
  plot.set_xaxis_label(xlabel);
  plot.set_yaxis_label(ylabel);
  
  return out;
}

template <typename W,
	  typename F,
	  typename A>
W FSEcorr(const W& w02m2Pi,const F& volDep,const A& arg,const double& power)
{
  return 1+w02m2Pi*volDep*exp(-arg)/pow(arg,power);
}

template <typename W,
	  typename F,
	  typename A>
W ansatz_Ldep(const W& w02m2Pi,const F& infVol,const F& volDep,const A& arg,const double& power)
{
  return infVol*FSEcorr(w02m2Pi,volDep,arg,power);
}

template <typename F>
void analysis_Ldep(const string& name,const vector<size_t>& ens40,const double& power,F getY)
{
  new_section("Studying "+name+" FSE");
  
  const size_t n40=ens40.size();
  
  vector<double> x(n40),plot_x(n40);
  djvec_t y(n40);
  jack_fit_t fit;
  djack_t infVol,volDep;
  const size_t iInfVol=fit.add_fit_par(infVol,"InfVol",0.141,0.001);
  const size_t iVolDep=fit.add_fit_par(volDep,"VolDep",3.96,0.001);
  
  for(size_t iens40=0;iens40<ens40.size();iens40++)
    {
      const size_t iens=ens40[iens40];
      const perens_t& ensUncorr=ensListUncorr[iens];
      perens_t& ensCorr=ensList[iens];
      const djack_t& amPion=ensUncorr.amPion;
      const size_t Lfra=ensUncorr.L;
      const djack_t w0_fr_a=ensUncorr.w0_fr_a;
      const djack_t arg=Lfra*amPion;
      const djack_t w02m2Pi=sqr((djack_t)(w0_fr_a*amPion));
      
      x[iens40]=exp(-arg.ave())/pow(arg.ave(),power);
      y[iens40]=getY(ensCorr);
      
      fit.add_point(y[iens40],[&power,w02m2Pi,arg,&iInfVol,&iVolDep](const vector<double>& p,const size_t ijack)
			      {
				return ansatz_Ldep(w02m2Pi[ijack],p[iInfVol],p[iVolDep],arg[ijack],power);
			      });
    }
  
  const auto status=fit.fit();
  // cout<<infVol.ave_err()<<endl;
  // cout<<volDep.ave_err()<<endl;
  
  const djack_t w02m2Pi_A40=sqr((djack_t)(ensListUncorr[ens40[0]].w0_fr_a*ensListUncorr[ens40[0]].amPion));
  grace_file_t plot("plots/ens40/"+name+".xmg");
  plot.write_vec_ave_err(x,y.ave_err());
  plot.set_transparency(0.3);
  plot.set_no_line();
  
  plot.write_polygon([&volDep,&infVol,&w02m2Pi_A40,power](const double& x)->djack_t
		     {
		       const djack_t arg=Brent_solve([y=x,power](double x){return exp(-x)/pow(x,power)-y;},1e-7,10);
		       
		       return ansatz_Ldep(w02m2Pi_A40,infVol,volDep,arg,power);
		     },1e-4,0.011,grace::RED);
  
  for(size_t iens=0;iens<ensList.size();iens++)
    {
      perens_t& ensUncorr=ensListUncorr[iens];
      perens_t& ensCorr=ensList[iens];
      const djack_t& amPion=ensUncorr.amPion;
      const size_t Lfra=ensUncorr.L;
      const djack_t w0_fr_a=ensUncorr.w0_fr_a;
      const djack_t arg=Lfra*amPion;
      const djack_t w02m2Pi=sqr((djack_t)(w0_fr_a*amPion));
      
      const djack_t corr=FSEcorr(w02m2Pi,volDep,arg,power),corrpercent=(corr-1)*100;
      const djack_t& yUncorr=getY(ensUncorr);
      djack_t& yCorr=getY(ensCorr);
      
      yCorr=yUncorr/corr;
      
      cout<<ensCorr.name<<"\t"<<yUncorr.ave_err()<<"\t"<<corrpercent.ave_err()<<" % \t"<<yCorr.ave_err()<<endl;
    }
  
  djvec_t yCorr(n40);
  for(size_t iens40=0;iens40<ens40.size();iens40++)
    {
      const size_t iens=ens40[iens40];
      perens_t& ensCorr=ensList[iens];
      yCorr[iens40]=getY(ensCorr);
    }
  plot.write_vec_ave_err(x,yCorr.ave_err());
  plot.set_no_line();
  
  plot.write_line([i=infVol.ave()](double x){return i;},1e-4,0.011);
}

int main()
{
  input_file_t analysis("analysis.txt");
  W0FrADir=analysis.read<string>("W0FrADir");
  afPionDir=analysis.read<string>("afPionDir");
  const size_t ext_njacks=analysis.read<size_t>("NJacks");
  set_njacks(ext_njacks);
  const size_t nEnsTot=analysis.read<size_t>("NEns");
  vector<size_t> nEnsPerBeta(nb,0);
  
  new_section("PREAMBLE");
  
  cout<<"List of exprimental quantities"<<endl;
  cout<<"MKBarExp: "<<MKaonBarExp<<" GeV"<<endl;
  cout<<"MEtaExp: "<<MEtaExp<<" GeV"<<endl;
  
  for(size_t iens=0;iens<nEnsTot;iens++)
    {
      /// Name of the ensemble
      const string name=analysis.read<string>();
      
      new_section("Working on ensemble: "+name);
      
      ensList.emplace_back(name);
      
      /// Last added ensemble
      perens_t& ens=ensList.back();
      
      // Add ml
      amlList[ens.ib].insert(ens.amq[0]);
      
      ens.
	loadW0fra().
	loadAfPion().
	loadCorrs().
	effectiveMasses().
	fitMasses().
	printMasses();
      
      nEnsPerBeta[ens.ib]++;
    }
  
  new_section("0.0040 analysis");
  
  // Count ensembles of 0.0040
  vector<size_t> ens40;
  for(size_t iens=0;iens<ensList.size();iens++)
    {
      const perens_t& ens=ensList[iens];
      if(ens.ib==0 and ens.amq[0]==0.0040)
	ens40.push_back(iens);
    }
  cout<<"N. of ensembles 0.0040: "<<ens40.size()<<endl;
  
  ensListUncorr=ensList;
  
  analysis_Ldep("w0FrA",ens40,1.0,[](perens_t& ens)->djack_t&{return ens.w0_fr_a;});
  analysis_Ldep("amPion",ens40,1.5,[](perens_t& ens)->djack_t&{return ens.amPion;});
  analysis_Ldep("amRho",ens40,1.5,[](perens_t& ens)->djack_t&{return ens.amPhi[0];});
  analysis_Ldep("afPion",ens40,1.5,[](perens_t& ens)->djack_t&{return ens.afPion;});
  analysis_Ldep("aPhiPion",ens40,1.5,[](perens_t& ens)->djack_t&{return ens.aPhiPion;});
  analysis_Ldep("amNucleon",ens40,1.0,[](perens_t& ens)->djack_t&{return ens.amNucleon;});
  for(size_t iq=1;iq<ensList[ens40.front()].nq;iq++)
    {
      analysis_Ldep("phi"+to_string(iq),ens40,1.5,[iq](perens_t& ens)->djack_t&{return ens.amPhi[iq];});
      analysis_Ldep("kaon"+to_string(iq),ens40,1.5,[iq](perens_t& ens)->djack_t&{return ens.amKaon[iq];});
      analysis_Ldep("omega"+to_string(iq),ens40,1.0,[iq](perens_t& ens)->djack_t&{return ens.amOmega[iq];});
    }
  
  new_section("ams for each ensemble");
  for(auto& ens : ensList)
    ens.
      setRatios().
      getAmsPhys();
  
  /// Function to interpolate
  auto interpolates=[](const string& name,const djvec_t data,const perens_t& ens,const djack_t& x)
		    {
		      return fitAndEval(ens.name+"/plots/"+name+".xmg",ens.amq,data,2,ens.amq[1]*0.9,ens.amq[3]*1.1,x);
		    };
  
  /// Interpolated ratioL to ams_phys for each ensemble
  djvec_t ratioLPerEns(nEnsTot);
  for(size_t iens=0;iens<nEnsTot;iens++)
    {
      const perens_t& ens=ensList[iens];
      
      ratioLPerEns[iens]=interpolates("ratioL",ens.ratioL,ens,ens.amsPhys);
    }
  
  ///////////////////// Determine ams beta per beta
  
  new_section("ams fit for each beta, with global slope");
  
  djack_t w0_alt;
  {
    jack_fit_t combo_fit;
    djvec_t combo_pars(3);
    combo_fit.add_fit_par(combo_pars[0],"w0M"+ssName,1,0.01);
    combo_fit.add_fit_par(combo_pars[1],"slope",0.1,0.01);
    combo_fit.add_fit_par(combo_pars[2],"slope_a2",0.0,0.01);
    
    auto F=[interpolates](const perens_t& ens) -> djack_t
	   {
	     return ens.w0_fr_a*interpolates(ssName,ens.amOnlyS,ens,ens.amsPhys);
	   };
    
    for(size_t iens=0;iens<nEnsTot;iens++)
      {
	const perens_t& ens=ensList[iens];
	combo_fit.add_point(F(ens),[&r=ratioLPerEns[iens],w0_fr_a=ens.w0_fr_a](const vector<double>& p,const size_t ijack)->double
				   {
				     return p[0]+r[ijack]*p[1]+p[2]/sqr(w0_fr_a[ijack]);
				   });
      }
    combo_fit.fit();
    cout<<"combo_fit_pars: "<<endl<<combo_pars.ave_err()<<endl;
    
    /// Plot determine phys for each beta
    grace_file_t plot("plots/w0FrAm"+ssName+"PhysAlone.xmg");
    vector<grace::color_t> colors{grace::TURQUOISE,grace::RED,grace::GREEN4,grace::ORANGE,grace::VIOLET,grace::BLACK};
    djvec_t physPerBeta(nb);
    for(size_t ib=0;ib<nb;ib++)
      {
	const perens_t* last_of_this_beta=nullptr;
	for(const perens_t& ens : ensList)
	  if(ens.ib==ib)
	    last_of_this_beta=&ens;
	if(last_of_this_beta==nullptr)
	  CRASH("No ensembleLoop found for beta %zu",ib);
	
	plot.write_polygon([&combo_pars,last_of_this_beta](const double& x)->djack_t
			   {
			     return combo_pars[0]+combo_pars[1]*x+combo_pars[2]/sqr(last_of_this_beta->w0_fr_a);
			   },0,0.060,colors[ib]);
	plot.new_data_set();
	plot.set_no_line();
	plot.set_all_colors(colors[ib]);
	plot.set_legend("\\xb\\0="+beta_tag[ib]);
	plot.set_no_line();
	
	for(size_t _iens=0;_iens<nEnsTot;_iens++)
	  {
	    const perens_t& ens=ensList[_iens];
	    
	    if(ens.ib==ib)
	      plot.write_ave_err(ratioLPerEns[_iens].ave(),F(ens).ave_err());
	  }
	
	{
	  physPerBeta[ib]=combo_pars[0]+ratioLExp*combo_pars[1]+combo_pars[2]/sqr(last_of_this_beta->w0_fr_a);
	}
	plot.write_ave_err(ratioLExp,physPerBeta[ib].ave_err());
	cout<<"(w0*mOnlyS)["<<ib<<"] "<<physPerBeta[ib].ave_err()<<endl;
	
	plot.write_polygon([&combo_pars](const double& x)->djack_t
			   {
			     return combo_pars[0]+combo_pars[1]*x;
			   },0,0.060,grace::VIOLET);
      }
    
    const djack_t w0_mOnlyS_exp=combo_pars[0]+ratioLExp;
    w0_alt=w0_mOnlyS_exp/MOnlySExp*0.197;
    cout<<"w0: "<<w0_alt<<" fm"<<endl;
  }
  
  /// Function to fit each quantity a
  auto fitQ=[nEnsTot,&ratioLPerEns](function<djack_t(const perens_t&)> F,const string& Qname,const bool WithCurvature)
	    {
	      jack_fit_t combo_fit;
	      djvec_t combo_pars(nb+3);
	      //djvec_t combo_pars(nb*2);
	      for(size_t i=0;i<nb;i++)
		combo_fit.add_fit_par(combo_pars[i],string("P")+to_string(i),F(ensList[0]).ave(),0.01);
	      combo_fit.add_fit_par(combo_pars[nb],"slope",0.1,0.01);
	      combo_fit.add_fit_par(combo_pars[nb+1],"slope_a2",0.0,0.01);
	      combo_fit.add_fit_par(combo_pars[nb+2],"curvature",0.0,0.01);
	      if(not WithCurvature)
		combo_fit.fix_par(nb+2);
		// {
		//   combo_fit.add_fit_par(combo_pars[2*i+0],string("P")+to_string(i),1,0.01);
		//   combo_fit.add_fit_par(combo_pars[2*i+1],string("Q")+to_string(i),1,0.01);
		// }
	      
	      for(size_t iens=0;iens<nEnsTot;iens++)
		{
		  const perens_t& ens=ensList[iens];
		  combo_fit.add_point(F(ens),[ib=ens.ib,&r=ratioLPerEns[iens],w0_fr_a=ens.w0_fr_a](const vector<double>& p,const size_t ijack)->double
					     {
					       const double& slope_a0=p[nb];
					       const double& slope_a2=p[nb+1];
					       const double& curvature=p[nb+2];
					       const double slope=slope_a0+slope_a2/sqr(w0_fr_a[ijack]);
					       
					       return p[ib]*(1+r[ijack]*slope+r[ijack]*r[ijack]*curvature);
					       //return p[2*ib+0]+p[2*ib+1]*r[ijack];
					     });
		}
	      combo_fit.fit();
	      cout<<"combo_fit_pars: "<<endl<<combo_pars.ave_err()<<endl;
	      
	      /// Plot determine phys for each beta
	      grace_file_t plot("plots/"+Qname+"_phys.xmg");
	      vector<grace::color_t> colors{grace::TURQUOISE,grace::RED,grace::GREEN4,grace::ORANGE,grace::VIOLET,grace::BLACK};
	      djvec_t physPerBeta(nb);
	      for(size_t ib=0;ib<nb;ib++)
		{
		  const perens_t* last_of_this_beta=nullptr;
		  for(const perens_t& ens : ensList)
		    if(ens.ib==ib)
		      last_of_this_beta=&ens;
		  if(last_of_this_beta==nullptr)
		    CRASH("No ensembleLoop found for beta %zu",ib);
		  
		  plot.write_polygon([&combo_pars,last_of_this_beta,ib](const double& x)->djack_t
		  //plot.write_polygon([&combo_pars,ib](const double& x)->djack_t
					 {
					   const djack_t& slope_a0=combo_pars[nb];
					   const djack_t& slope_a2=combo_pars[nb+1];
					   const djack_t& curvature=combo_pars[nb+2];
					   const djack_t slope=slope_a0+slope_a2/sqr(last_of_this_beta->w0_fr_a);
					   
					   return combo_pars[ib]*(1+x*slope+x*x*curvature);
					   //return combo_pars[2*ib+0]+combo_pars[2*ib+1]*x;
					 },0,0.160,colors[ib]);
		  plot.new_data_set();
		  plot.set_no_line();
		  plot.set_all_colors(colors[ib]);
		  plot.set_legend("\\xb\\0="+beta_tag[ib]);
		  plot.set_no_line();
		  
		  for(size_t _iens=0;_iens<nEnsTot;_iens++)
		    {
		      const perens_t& ens=ensList[_iens];
		      
		      if(ens.ib==ib)
			plot.write_ave_err(ratioLPerEns[_iens].ave(),F(ens).ave_err());
		    }
		  
		  physPerBeta[ib]=combo_pars[ib]*(1+ratioLExp*combo_pars[nb]+ratioLExp*ratioLExp*combo_pars[nb+2]);
		  //physPerBeta[ib]=combo_pars[2*ib+0]+ratioLExp*combo_pars[2*ib+1];
		  plot.write_ave_err(ratioLExp,physPerBeta[ib].ave_err());
		  cout<<Qname<<"["<<ib<<"] "<<physPerBeta[ib].ave_err()<<endl;
		}
	      
	      return physPerBeta;
	    };
  
  const bool WithCurvature=true,NoCurvature=false;
  const djvec_t aPerBeta=fitQ([w0_alt](const perens_t& ens)->djack_t{return w0_alt/ens.w0_fr_a;},"a",NoCurvature);
  const djvec_t amsPhysPerBeta=fitQ([](const perens_t& ens)->djack_t{return ens.amsPhys;},"ams",NoCurvature);
  // const djvec_t amsPhiPhysPerBeta=fitQ([](const perens_t& ens)->djack_t{return ens.amsPhiPhys;},"amsPhi",NoCurvature);
  const djvec_t amlFrRatioLPhysPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return ens.amq[0]/interpolates("ratioL",ens.ratioL,ens,ens.amsPhys);},"aml",WithCurvature);
  const djvec_t W0FrA_amOnlySPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return ens.w0_fr_a*interpolates(ssName,ens.amOnlyS,ens,ens.amsPhys);},"w0_fr_a_amOnlyS",NoCurvature);
  const djvec_t amOnlySPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return interpolates(ssName,ens.amOnlyS,ens,ens.amsPhys);},ssName,NoCurvature);
  const djvec_t amOmegaFrAmPhiPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return interpolates("amOmegaFrAmPhi",ens.amOmegaFrAmPhi,ens,ens.amsPhys);},"amOmegaFrAmPhi",NoCurvature);
  const djvec_t amNucleonFrAmOnlySPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return ens.amNucleon/interpolates(ssName,ens.amOnlyS,ens,ens.amsPhys);},"amNucleonFrAm"+ssName,NoCurvature);
  const djvec_t afPionFrAmOnlySPerBeta=fitQ([interpolates](const perens_t& ens)->djack_t{return ens.afPion/interpolates(ssName,ens.amOnlyS,ens,ens.amsPhys);},"afPionFrAm"+ssName,WithCurvature);
  const djvec_t amlPhysPerBeta=amlFrRatioLPhysPerBeta*ratioLExp;
  const djvec_t a=amOnlySPerBeta/MOnlySExp,ainfm=a*0.197;
  cout<<"a: "<<endl<<ainfm.ave_err()<<" fm"<<endl;
  
  const djvec_t amsFrAmlPhysPerBeta=amsPhysPerBeta/amlPhysPerBeta;
  cout<<"ams/aml: "<<endl<<amsFrAmlPhysPerBeta.ave_err()<<endl;
  // const djvec_t amsPhysFrPhysPhiPerBeta=amsPhysPerBeta/amsPhiPhysPerBeta;
  // cout<<"ams/amsPhi: "<<endl<<amsPhysFrPhysPhiPerBeta.ave_err()<<endl;
  
  djvec_t a_lit(3);
  a_lit[0].fill_gauss(0.0886/0.197,0.0027/0.197,2142);
  a_lit[1].fill_gauss(0.0815/0.197,0.0021/0.197,2143);
  a_lit[2].fill_gauss(0.0619/0.197,0.0011/0.197,2144);
  cout<<"a lit: "<<a_lit.ave_err()<<endl;
  
  const vector<vector<ave_err_t>> Zp_ae({{{0.529,0.007},{0.509,0.004},{0.516,0.002}},{{0.574,0.004},{0.546,0.002},{0.545,0.002}}});
  // const double conva=0.7996,convb=0.8330,convc=0.9263;
  // const vector<vector<ave_err_t>> Zp_ae({{{0.421613/conva,0.004362/conva},{0.437983/convb,0.004253/convb},{0.498469/convc,0.001839/convc}},
  // 					 {{0.446678/conva,0.002027/conva},{0.459713/convb,0.002155/convb},{0.513690/convc,0.000997/convc}}});
  const vector<vector<ave_err_t>> Za_ae({{{0.731,0.008},{0.737,0.005},{0.762,0.004}},{{0.703,0.002},{0.714,0.002},{0.752,0.002}}});
  const int Zp_seed[nb]={82223,224335,34434};
  djvec_t Zp(nb);
  
  for(size_t ibeta=0;ibeta<nb;ibeta++)
    Zp[ibeta].fill_gauss(Zp_ae[1][ibeta],Zp_seed[ibeta]);
  
  //double Zp[]={0.529,0.509,0.516};
  vector<double> a2(nb);
  djvec_t msRenPerBeta(nb);
  for(size_t ib=0;ib<nb;ib++)
    {
      msRenPerBeta[ib]=amsPhysPerBeta[ib]/a[ib]/Zp[ib];
      a2[ib]=sqr(a[ib]).ave();
    }
  cout<<"ms: "<<endl<<msRenPerBeta.ave_err()<<endl;
  
  const djvec_t msRenPars=poly_fit(a2,msRenPerBeta,1,0,0.3,"plots/msRen.xmg");
  const djack_t msRen=msRenPars[0];
  cout<<"ms: "<<smart_print(msRen)<<endl;
  
  const djvec_t amNucleonFrAmOnlySPars=poly_fit(a2,amNucleonFrAmOnlySPerBeta,1,0,0.3,"plots/mNucleonFrAm"+ssName+".xmg");
  const djack_t mNucleon=amNucleonFrAmOnlySPars[0]*MOnlySExp;
  cout<<"mNucleon: "<<smart_print(mNucleon)<<" GeV, exp: "<<MProtonExp<<" GeV"<<endl;
  
  const djvec_t amOmegaFrAmPhiPars=poly_fit(a2,amOmegaFrAmPhiPerBeta,1,0,0.3,"plots/amOmegaFrAmPhi.xmg");
  const djack_t mOmegaFrMPhi=amOmegaFrAmPhiPars[0];
  cout<<"mOmega/mPhi: "<<smart_print(mOmegaFrMPhi)<<" , exp: "<<MOmegaFrMPhiExp<<endl;
  
  const djvec_t W0FrA_amOnlyS=poly_fit(a2,W0FrA_amOnlySPerBeta,1,0,0.3,"plots/W0FrA_am"+ssName+".xmg");
  const djack_t W0=W0FrA_amOnlyS[0]/MOnlySExp*0.197;
  cout<<"W0: "<<smart_print(W0)<<" fm, exp: "<<W0Exp<<" fm"<<endl;
  
  const djvec_t afPion_fr_amOnlyS=poly_fit(a2,afPionFrAmOnlySPerBeta,1,0,0.3,"plots/afPionFrAm"+ssName+".xmg");
  const djack_t fPion=afPion_fr_amOnlyS[0]*MOnlySExp;
  cout<<"fPi: "<<smart_print(fPion)<<" GeV, exp: "<<fPionExp<<" fm"<<endl;
  
  const djvec_t amsFrAmlPars=poly_fit(a2,amsFrAmlPhysPerBeta,1,0,0.3,"plots/amsFrAml.xmg");
  const djack_t amsFrAml=amsFrAmlPars[0];
  cout<<"ams/aml: "<<smart_print(amsFrAml)<<", exp: "<<msFrMlExp<<endl;
  
  const djack_t mlRen=msRen/amsFrAml;
  cout<<"ml: "<<smart_print(mlRen)<<endl;
  
  return 0;
}
