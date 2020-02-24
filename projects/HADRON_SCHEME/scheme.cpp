#include <tranalisi.hpp>

#include <set>

const size_t nb=3;

/// List of ml ordered by ib
map<size_t,set<double>> amlList;

const double MPionZeroExp=0.13498;
const double MKaonPlusExp=0.49368;
const double MKaonZeroExp=0.49761;
const double MKaonBarExp=sqrt((sqr(MKaonPlusExp)+sqr(MKaonZeroExp))/2);
const double MEtaExp=sqrt(2*sqr(MKaonBarExp)-sqr(MPionZeroExp));
const double MOmegaExp=1.672;
const double MProtonExp=0.938;
const double ratioAExp=MOmegaExp/MEtaExp;
const double ratioBExp=sqr(MPionZeroExp/MOmegaExp);

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
  
  /// Reads a meson
  djvec_t load_meson(const size_t iq1,const size_t iq2) const
  {
    return read_djvec(combine("%s/jacks/mes_contr_P5P5_SM_q%zu__SM_q%zu",name.c_str(),iq1,iq2),T).symmetrized();
  }
  
  djvec_t cPion,ePion;
  djack_t amPion;
  
  djvec_t cNucleon,eNucleon;
  djack_t amNucleon;
  
  vector<djvec_t> cKaon,eKaon;
  djvec_t amKaon;
  
  vector<djvec_t> cOmega,eOmega;
  djvec_t amOmega;
  
  djvec_t amEtaSS;
  
  djvec_t ratioA;
  
  djvec_t ratioB;
  
  djack_t amsPhys;
  
  /// Load the correlators
  perens_t& loadCorrs()
  {
    cPion=load_meson(0,0);
    cNucleon=load_nucleon(0);
    
    for(auto& c : {&cKaon,&cOmega})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	cOmega[iq]=load_omega(iq);
	cKaon[iq]=load_meson(0,iq);
      }
    
    return *this;
  }
  
  /// Computes the effective masses
  perens_t& effectiveMasses()
  {
    ePion=effective_mass(cPion);
    eNucleon=effective_mass(cNucleon);
    
    for(auto& c : {&eKaon,&eOmega})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	eOmega[iq]=effective_mass(cOmega[iq],TH,-1);
	eKaon[iq]=effective_mass(cKaon[iq]);
	
      }
    
    return *this;
  }
  
  /// Fit the effective masses
  perens_t& fitMasses()
  {
    amPion=constant_fit(ePion,tminPion,tmaxPion,name+"/plots/pion.xmg");
    amNucleon=constant_fit(eNucleon,tminNucleon,tmaxNucleon,name+"/plots/nucleon.xmg");
    
    for(auto& c : {&amKaon,&amOmega,&amEtaSS,&ratioA,&ratioB})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	amOmega[iq]=constant_fit(eOmega[iq],tminOmega,tmaxOmega,name+combine("/plots/omega_%zu.xmg",iq));
	amKaon[iq]=constant_fit(eKaon[iq],tminKaon,tmaxKaon,name+combine("/plots/kaon_%zu.xmg",iq));
	amEtaSS[iq]=sqrt(2*sqr(amKaon[iq])-sqr(amPion));
	ratioA[iq]=amOmega[iq]/amEtaSS[iq];
	ratioB[iq]=sqr((djack_t(amPion/amOmega[iq])));
      }
    
    return *this;
  }
  
  /// Print fitted masses
  void printMasses()
  {
    cout<<"aMNuc: "<<amNucleon.ave_err()<<endl;
    for(size_t iq=0;iq<nq;iq++)
      {
	cout<<"MK: "<<amKaon[iq].ave_err()<<endl;
	cout<<"MEtaSS: "<<amEtaSS[iq].ave_err()<<endl;
	cout<<"MOmega: "<<amOmega[iq].ave_err()<<endl;
      }
  }
  
  /// Computes ams phys
  perens_t& getAmsPhys()
  {
    grace_file_t ratioAPlot(name+"/plots/ratioA.xmg");
    
    for(size_t iq=1;iq<nq;iq++)
      ratioAPlot.write_ave_err(amq[iq],ratioA[iq].ave_err());
    
    /// Fit the ratioA with a second order polynomial
    const djvec_t ratioAPars=poly_fit(amq,ratioA,2,amq[1]*0.9);
    ratioAPlot.write_polygon([ratioAPars](const double x){return poly_eval(ratioAPars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::GREEN4);
    const djack_t amsQuad=parab_solve(ratioAPars,ratioAExp);
    ratioAPlot.write_ave_err(amsQuad.ave_err(),{ratioAExp,0.0});
    cout<<"ams(quad): "<<amsQuad.ave_err()<<endl;
    
    /// Fit the ratioA at all ms with a first order polynomial
    const djvec_t ratioALin3Pars=poly_fit(amq,ratioA,1,amq[1]*0.9);
    ratioAPlot.write_polygon([ratioALin3Pars](const double x){return poly_eval(ratioALin3Pars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    const djack_t amsLin3=-(ratioALin3Pars[0]-ratioAExp)/ratioALin3Pars[1];
    cout<<"ams(lin3): "<<amsLin3.ave_err()<<endl;
    
    /// Fit the ratioA at all ms with a first order polynomial
    const djvec_t ratioALin2Pars=poly_fit(amq,ratioA,1,amq[2]*0.9);
    ratioAPlot.write_polygon([ratioALin2Pars](const double x){return poly_eval(ratioALin2Pars,x);},amq[2]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    const djack_t amsLin2=-(ratioALin2Pars[0]-ratioAExp)/ratioALin2Pars[1];
    cout<<"ams(lin2): "<<amsLin2.ave_err()<<endl;
    
    amsPhys=amsQuad;
    
    return *this;
  }
};

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
djack_t fitAndEval(const string& path,const vector<double>& x,const djvec_t& y,const size_t order,const double min,const double max,const T x0)
{
  const djvec_t pars=poly_fit(x,y,order,min,max);
  const djack_t out=poly_eval(pars,x0);
  
  grace_file_t plot(path);
  plot.write_polygon([&pars](const double& x){return poly_eval(pars,x);},min,max);
  plot.set_no_line();
  plot.write_vec_ave_err(x,y.ave_err());
  plot.set_all_colors(grace::GREEN4);
  plot.write_ave_err(fitEvalX(x0),out.ave_err());
  
  return out;
}

int main()
{
  input_file_t analysis("analysis.txt");
  const size_t ext_njacks=analysis.read<size_t>("NJacks");
  set_njacks(ext_njacks);
  const size_t nensTot=analysis.read<size_t>("NEns");
  vector<size_t> nensPerBeta(nb,0);
  
  cout<<"MKBar: "<<MKaonBarExp<<" GeV"<<endl;
  cout<<"MEtaExp: "<<MEtaExp<<" GeV"<<endl;
  
  vector<perens_t> ensList;
  
  for(size_t iens=0;iens<nensTot;iens++)
    {
      const string name=analysis.read<string>();
      cout<<name<<endl;
      
      ensList.emplace_back(name);
      
      /// Gets last added ensemble
      perens_t& ens=ensList.back();
      
      // Add ml
      amlList[ens.ib].insert(ens.amq[0]);
      
      ens.
	loadCorrs().
	effectiveMasses().
	fitMasses().
	getAmsPhys();
      
      nensPerBeta[ens.ib]++;
    }
  
  // Analysis beta per beta
  djvec_t a(nb);
  djvec_t ratioBPerEns(nensTot);
  for(size_t ib=0;ib<nb;ib++)
    {
      /// Function to interpolate
      auto interpolates=[](const string& name,const djvec_t data,const perens_t& ens)
			{
			  return fitAndEval(ens.name+"/plots/"+name+".xmg",ens.amq,data,2,ens.amq[1]*0.9,ens.amq[3]*1.1,ens.amsPhys);
			};
      
      const size_t& nens=nensPerBeta[ib];
      
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<"Beta "<<ib<<endl;
      
      if(nens<=1)
	cout<<" Skipping"<<endl;
      else
	{
	  /// Slice ensembles
	  vector<double> amlPerEns(nens);
	  djvec_t amPionPerEns(nens);
	  djvec_t amsPhysPerEns(nens),ratioBPerEnsSlice(nens),amOmegaPerEns(nens);
	  size_t iens=0;
	  for(size_t _iens=0;_iens<nensTot;_iens++)
	    {
	      const perens_t& ens=ensList[_iens];
	      
	      if(ens.ib==ib)
		{
		  amlPerEns[iens]=ens.amq[0];
		  amPionPerEns[iens]=ens.amPion;
		  amsPhysPerEns[iens]=ens.amsPhys;
		  
		  /// Interpolates ratioB and amOmega at phys ams
		  ratioBPerEns[_iens]=ratioBPerEnsSlice[iens]=interpolates("ratioB",ens.ratioB,ens);
		  amOmegaPerEns[iens]=interpolates("omega",ens.amOmega,ens);
		  
		  iens++;
		}
	    }
	  
	  auto extrapolates=[&ib,&ratioBPerEnsSlice](const string& name,const djvec_t data)
			    {
			      return fitAndEval("plots/"+name+"_beta"+to_string(ib)+".xmg",ratioBPerEnsSlice.ave(),data,1,0.0,0.065,ratioBExp);
			    };
	  
	  /// Determine the strange mass
	  const djack_t amsPhys=extrapolates("amsPhys",amsPhysPerEns);
	  cout<<"ams: "<<amsPhys.ave_err()<<endl;
	  
	  /// Determine the omega mass
	  const djack_t amOmega=extrapolates("amOmega",amOmegaPerEns);
	  cout<<"amOmega: "<<amOmega.ave_err()<<endl;
	  
	  /// Determine the lattice spacing
	  a[ib]=amOmega/MOmegaExp;
	  const djack_t aInv=1/a[ib];
	  cout<<"a^-1: "<<smart_print(aInv)<<" GeV^-1"<<endl;
	}
    }
  
  plan_fit_data_t<djack_t> contChirFitData;
  djvec_t mNucleonPerEns(nensTot);
  for(size_t iens=0;iens<nensTot;iens++)
    {
      const perens_t& ens=ensList[iens];
      const size_t ib=ens.ib;
      
      vector<double> x(3);
      x[0]=1.0;
      x[1]=ratioBPerEns[iens].ave();
      x[2]=sqr(a[ens.ib].ave());
      
      mNucleonPerEns[iens]=ens.amNucleon/a[ib];
      contChirFitData.push_back(make_tuple(x,mNucleonPerEns[iens]));
    }
  
  const auto coeffs=plan_fit(contChirFitData);
  const vector<double> contChirPoint{1.0,ratioBExp,0.0};
  const djack_t MN=plan_eval(coeffs,contChirPoint);
  
  cout<<"MNucleon: "<<smart_print(MN)<<" GeV"<<endl;
  
  grace_file_t plotMNucleonFit("plots/MN.xmg");
  vector<grace::color_t> colors{grace::TURQUOISE,grace::RED,grace::GREEN4,grace::ORANGE,grace::VIOLET,grace::BLACK};
  for(size_t ib=0;ib<nb;ib++)
    if(nensPerBeta[ib]>1)
      {
	plotMNucleonFit.write_polygon([&coeffs,&a,&ib](const double& x){return plan_eval(coeffs,vector<double>{1.0,x,sqr(a[ib].ave())});},0.0,0.065,colors[ib]);
	
	plotMNucleonFit.set_all_colors(colors[ib]);
	plotMNucleonFit.set_no_line();
	for(size_t _iens=0;_iens<nensTot;_iens++)
	  if(ensList[_iens].ib==ib)
	    plotMNucleonFit.write_ave_err(ratioBPerEns[_iens].ave(),mNucleonPerEns[_iens].ave_err());
	plotMNucleonFit.new_data_set();
      }
  plotMNucleonFit.write_ave_err(ratioBExp,MN.ave_err());
  plotMNucleonFit.set_all_colors(colors[nb]);
  plotMNucleonFit.write_polygon([&coeffs](const double& x){return plan_eval(coeffs,vector<double>{1.0,x,0.0});},0.0,0.065,colors[nb]);

  plotMNucleonFit.write_ave_err(ratioBExp,{MProtonExp,0.0});
  plotMNucleonFit.set_all_colors(colors[nb+1]);
  
  // /// Compute MN, MOmega and M^2Pi/ml
  // djvec_t amlFrMsPerEns(nens);
  // djvec_t MNPerEns(nens);
  // djvec_t M2PiPerEns(nens);
  // djvec_t M2PiFrAmlAmsPerEns(nens);
  // djvec_t MOmega2Perens(nens);
  // djvec_t RatioA2PerEns(nens);
  // for(size_t iens=0;iens<nens;iens++)
  //   {
  //     perens_t& ens=ensList[iens];
      
  //     const double aml=ens.amq[0];
  //     const djack_t MPi=ens.amPion*aInv;
  //     amlFrMsPerEns[iens]=aml/amsPhys;
  //     MNPerEns[iens]=ens.amNucleon*aInv;
  //     M2PiFrAmlAmsPerEns[iens]=sqr(MPi)*amsPhys/aml;
  //     M2PiPerEns[iens]=sqr(MPi);
  //     MOmega2Perens[iens]=ens.amOmega[2]*aInv;
  //     RatioA2PerEns[iens]=ens.ratioA[2];
  //   }
  
  // /// Fit M^2Pi/ml as a function of aml
  // grace_file_t plotM2PiFrAml("plots/M2PiFrAmlAms.xmg");
  // plotM2PiFrAml.write_vec_ave_err(amlFrMsPerEns.ave(),M2PiFrAmlAmsPerEns.ave_err());
  // const djvec_t M2PiFrAmlAmsFitPars=poly_fit(amlFrMsPerEns.ave(),M2PiFrAmlAmsPerEns,1);
  // plotM2PiFrAml.write_line([](const double x){return sqr(MPionZeroExp)/x;},0.03,0.06,grace::GREEN4);
  // plotM2PiFrAml.write_polygon([&M2PiFrAmlAmsFitPars](const double x) -> djack_t{return poly_eval(M2PiFrAmlAmsFitPars,x);},0,0.35,grace::VIOLET);
  
  // /// Plots M^2Pi/ml as a function of exp(-MPi*L)
  // grace_file_t plotAmPionFunVol("plots/amPionFunvol.xmg");
  // grace_file_t plotM2PiFrAmlFunVol("plots/M2PiFrAmlAmsFunVol.xmg");
  // grace_file_t plotAmsPhysFunVol("plots/amsPhysFunVol.xmg");
  // grace_file_t plotMNFunVol("plots/MNFunVol.xmg");
  // grace_file_t plotMOmega2FunVol("plots/MOmega2FunVol.xmg");
  // grace_file_t plotRatio2FunVol("plots/Ratio2FunVol.xmg");
  // vector<tuple<grace_file_t*,djvec_t*>> plotsFunVol={make_tuple(&plotM2PiFrAmlFunVol,&M2PiFrAmlAmsPerEns),
  // 						     make_tuple(&plotAmsPhysFunVol,&amsPhysPerEns),
  // 						     make_tuple(&plotAmPionFunVol,&amPionPerEns),
  // 						     make_tuple(&plotMNFunVol,&MNPerEns),
  // 						     make_tuple(&plotMOmega2FunVol,&MOmega2Perens),
  // 						     make_tuple(&plotRatio2FunVol,&RatioA2PerEns)};
  // for(auto& i : plotsFunVol)
  //   {
  //     grace_file_t& p=*get<grace_file_t*>(i);
  //     p.set_all_colors(grace::MAROON);
  //     p.set_color_scheme({grace::TURQUOISE,grace::RED,grace::GREEN4,grace::ORANGE,grace::VIOLET,grace::BLACK});
  //   }
  
  // for(size_t ib=0;ib<nb;ib++)
  //   for(const double& aml : amlList[ib])
  //     {
  // 	// Loop over all ensembles to gather all with given ib and aml
  // 	for(size_t iens=0;iens<nens;iens++)
  // 	  {
  // 	    const perens_t& ens=ensList[iens];
	    
  // 	    if(ens.ib==ib and ens.amq[0]==aml)
  // 	      {
  // 		const djack_t X=exp(-ens.amPion*ens.L);
  // 		for(auto& i : plotsFunVol)
  // 		  {
  // 		    grace_file_t& p=*get<grace_file_t*>(i);
  // 		    const djvec_t& y=*get<djvec_t*>(i);
  // 		    p<<"#"<<ens.name<<"\n";
  // 		    p.write_ave_err(X.ave_err(),y[iens].ave_err());
  // 		  }
  // 	      }
  // 	  }
	
  // 	for(auto& i : plotsFunVol)
  // 	  {
  // 	    grace_file_t& p=*get<grace_file_t*>(i);
  // 	    p.set_legend(combine("%zu_%lg",ib,aml));
  // 	    p.new_data_set();
  // 	  }
  //     }
  
  // /// Solve for aml
  // djvec_t M2PiFitPars(3);
  // M2PiFitPars[0]=0;
  // M2PiFitPars[1]=M2PiFrAmlFitPars[0];
  // M2PiFitPars[2]=M2PiFrAmlFitPars[1];
  // const djack_t amlPhys=parab_solve(M2PiFitPars,sqr(MPionZeroExp),false);
  // cout<<"aml: "<<smart_print(amlPhys)<<endl;
  
  // /// Plots M2pi
  // grace_file_t plotM2Pi("plots/M2Pi.xmg");
  // plotM2Pi.write_polygon([&M2PiFitPars](const double x) -> djack_t{return poly_eval(M2PiFitPars,x);},0,0.0110,grace::VIOLET);
  // plotM2Pi.write_vec_ave_err(aml,M2Pi.ave_err());
  // plotM2Pi.write_line([](const double x){return sqr(MPionZeroExp);},0.0007,0.0014,grace::GREEN4);
  // plotM2Pi.write_ave_err(amlPhys.ave_err(),{sqr(MPionZeroExp),0});
  
  // grace_file_t plotMN("plots/MN.xmg");
  // djvec_t MNFitPars=poly_fit(M2Pi.ave(),MN,1);;
  // plotMN.write_polygon([&MNFitPars](const double x) -> djack_t{return poly_eval(MNFitPars,x);},0.0,0.250,grace::VIOLET);
  // //plotMN.write_line([](const double x){return MProton;},0.0,0.0110,grace::GREEN4);
  // const djack_t MNextr=poly_eval(MNFitPars,sqr(MPionZeroExp));
  // plotMN.write_vec_ave_err(M2Pi.ave(),MN.ave_err());
  // plotMN.write_ave_err(sqr(MPionZeroExp),MNextr.ave_err());
  // cout<<"MN: "<<smart_print(MNextr)<<" GeV (phys: "<<MProtonExp<<" GeV)"<<endl;
  // /// Finds the physical ams
    
    // const djack_t lata=1/lata_inv;
    // cout<<"latinv "<<lata_inv.ave_err()<<" GeV"<<endl;
    // cout<<"lat "<<lata.ave_err()<<" GeV^-1"<<endl;
    
    // const double rat2_phys=MPionZero/MOmega;
    // const djack_t rat2=mpion/momega_inte/rat2_phys;
    // cout<<"rat2: "<<sqr(rat2).ave_err()<<endl;
  return 0;
}
