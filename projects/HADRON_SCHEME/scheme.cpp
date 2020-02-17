#include <tranalisi.hpp>

const double MPionZero=0.13498;
const double MKaonPlus=0.49368;
const double MKaonZero=0.49761;
const double MKaonBar=sqrt((sqr(MKaonPlus)+sqr(MKaonZero))/2);
const double MEtaExp=sqrt(2*sqr(MKaonBar)-sqr(MPionZero));
const double MOmega=1.672;
const double ratioExp=MOmega/MEtaExp;

/// Solve for the x of a parabola passing through d
djack_t parab_solve(const djvec_t& pars,const djack_t d)
{
  const djack_t& a=pars[2],&b=pars[1],&c=pars[0]-ratioExp;
  const djack_t& x0=(-b-sqrt(b*b-4*a*c))/(2*a);
  
  return x0;
}

struct perens_t
{
  string name;
  size_t nq;
  vector<double> amq;
  size_t T,TH;
  size_t tmin,tmax;
  
  perens_t(const string& name) : name(name)
  {
    input_file_t input(name+"/input");
    
    T=input.read<size_t>("T");
    TH=T/2;
    nq=input.read<int>("amq");
    amq.resize(nq);
    for(size_t iq=0;iq<nq;iq++)
      amq[iq]=input.read<double>();
    
    tmin=input.read<size_t>("tmin");
    tmax=input.read<size_t>("tmax");
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
  
  vector<djvec_t> eEtaSS;
  djvec_t amEtaSS;
  
  vector<djvec_t> cRatio;
  djvec_t ratio;
  
  djack_t amsLin2;
  djack_t amsLin3;
  djack_t amsQuad;
  
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
    
    for(auto& c : {&eKaon,&eOmega,&eEtaSS,&cRatio})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	eOmega[iq]=effective_mass(cOmega[iq]);
	eKaon[iq]=effective_mass(cKaon[iq]);
	eEtaSS[iq]=sqrt(2*sqr(eKaon[iq])-sqr(ePion[iq]));
	cRatio[iq]=eOmega[iq]/eEtaSS[iq];
      }
    
    return *this;
  }
  
  /// Fit the effective masses
  perens_t& fitMasses()
  {
    amPion=constant_fit(ePion,tmin,tmax,name+"/plots/pion.xmg");
    amNucleon=constant_fit(eNucleon,tmin,tmax,name+"/plots/nucleon.xmg");
    
    for(auto& c : {&amKaon,&amOmega,&amEtaSS,&ratio})
      c->resize(nq);
    
    for(size_t iq=0;iq<nq;iq++)
      {
	amOmega[iq]=constant_fit(eOmega[iq],tmin,tmax,name+combine("/plots/omega_%zu.xmg",iq));
	amKaon[iq]=constant_fit(eKaon[iq],tmin,tmax,name+combine("/plots/kaon_%zu.xmg",iq));
	amEtaSS[iq]=constant_fit(eEtaSS[iq],tmin,TH-3,name+combine("/plots/etaSS_%zu.xmg",iq));
	ratio[iq]=constant_fit(cRatio[iq],tmin,tmax,name+combine("/plots/ratio_%zu.xmg",iq));
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
    grace_file_t ratioPlot(name+"/plots/ratio.xmg");
    
    for(size_t iq=1;iq<nq;iq++)
      ratioPlot.write_ave_err(amq[iq],ratio[iq].ave_err());
    
    /// Fit the ratio with a second order polynomial
    const djvec_t ratioPars=poly_fit(amq,ratio,2,amq[1]*0.9);
    ratioPlot.write_polygon([ratioPars](const double x){return poly_eval(ratioPars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::GREEN4);
    amsQuad=parab_solve(ratioPars,ratioExp);
    ratioPlot.write_ave_err(amsQuad.ave_err(),{ratioExp,0.0});
    cout<<"ams(quad): "<<amsQuad.ave_err()<<endl;
    
    /// Fit the ratio at all ms with a first order polynomial
    const djvec_t ratioLin3Pars=poly_fit(amq,ratio,1,amq[1]*0.9);
    ratioPlot.write_polygon([ratioLin3Pars](const double x){return poly_eval(ratioLin3Pars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    amsLin3=-(ratioLin3Pars[0]-ratioExp)/ratioLin3Pars[1];
    cout<<"ams(lin3): "<<amsLin3.ave_err()<<endl;
    
    /// Fit the ratio at all ms with a first order polynomial
    const djvec_t ratioLin2Pars=poly_fit(amq,ratio,1,amq[2]*0.9);
    ratioPlot.write_polygon([ratioLin2Pars](const double x){return poly_eval(ratioLin2Pars,x);},amq[2]*0.9,amq[nq-1]*1.1,grace::VIOLET);
    amsLin2=-(ratioLin2Pars[0]-ratioExp)/ratioLin2Pars[1];
    cout<<"ams(lin2): "<<amsLin2.ave_err()<<endl;
    
    return *this;
  }
  
    // /// Finds the physical ams
    // ratio_plot.set_all_colors(grace::ORANGE);
    // ratio_plot.write_ave_err(ams_phys.ave_err(),{ratioExp,0.0});
    
    // const djack_t momega_inte=poly_eval(momegaPars,ams_phys);
    
    // const djack_t lata_inv=MOmega/momega_inte;
    // const djack_t lata=1/lata_inv;
    // cout<<"latinv "<<lata_inv.ave_err()<<" GeV"<<endl;
    // cout<<"lat "<<lata.ave_err()<<" GeV^-1"<<endl;
    
    // const double rat2_phys=MPionZero/MOmega;
    // const djack_t rat2=mpion/momega_inte/rat2_phys;
    // cout<<"rat2: "<<sqr(rat2).ave_err()<<endl;
};

int main()
{
  input_file_t analysis("analysis.txt");
  const size_t ext_njacks=analysis.read<size_t>("NJacks");
  set_njacks(ext_njacks);
  const size_t nens=analysis.read<size_t>("NEns");
  
  cout<<"MKBar: "<<MKaonBar<<" GeV"<<endl;
  cout<<"MEtaExp: "<<MEtaExp<<" GeV"<<endl;
  cout<<"MEtaExp: "<<ratioExp<<endl;
  
  vector<perens_t> ensList;

  for(size_t iens=0;iens<nens;iens++)
    {
      const string name=analysis.read<string>();
      ensList.emplace_back(name);
      
      cout<<name<<endl;
      
      ensList.back().
	loadCorrs().
	effectiveMasses().
	fitMasses().
	getAmsPhys();
    }
  
  grace_file_t plotAms("plots/ams.xmg");
  for(auto& ens : ensList)
    {
      plotAms.write_ave_err(ens.amq[0],ens.amsLin3.ave_err());
    }
  
  return 0;
}
