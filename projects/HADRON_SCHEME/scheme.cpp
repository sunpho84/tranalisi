#include <tranalisi.hpp>

size_t nq;
vector<double> amq;
size_t T,TH;

const double MPionZero=0.13498;
const double MKaonPlus=0.49368;
const double MKaonZero=0.49761;
const double MKaonBar=sqrt((sqr(MKaonPlus)+sqr(MKaonZero))/2);
const double MEtaExp=sqrt(2*sqr(MKaonBar)-sqr(MPionZero));
const double MOmega=1.672;
const double ratioExp=MOmega/MEtaExp;

size_t tmin,tmax;

djvec_t load_bar(const string& tag,const int iproj,const int iWick)
{
  string path=combine("jacks/bar_alt_contr_%s_proj_%d_Wick_%d",tag.c_str(),iproj,iWick);
  
  return -read_djvec(path,T).symmetrized(-1);
}

string bar_tag(const size_t iq)
{
  return combine("SM_q%zu_q%zu_q%zu",iq,iq,iq);
}

djvec_t load_nucleon(const size_t iq)
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

djvec_t load_omega(const size_t iq)
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

djvec_t load_meson(const size_t iq1,const size_t iq2)
{
  return read_djvec(combine("jacks/mes_contr_P5P5_SM_q%zu__SM_q%zu",iq1,iq2),T).symmetrized();
}

djack_t parab_solve(const djvec_t& pars,const djack_t d)
{
  const djack_t& a=pars[2],&b=pars[1],&c=pars[0]-ratioExp;
  const djack_t& x0=(-b-sqrt(b*b-4*a*c))/(2*a);
  
  return x0;
}

int main()
{
  input_file_t input("input");
  T=input.read<size_t>("T");
  TH=T/2;
  nq=input.read<int>("amq");
  amq.resize(nq);
  for(size_t iq=0;iq<nq;iq++)
    amq[iq]=input.read<double>();
  size_t ext_njacks=input.read<size_t>("NJacks");
  set_njacks(ext_njacks);
  tmin=input.read<size_t>("tmin");
  tmax=input.read<size_t>("tmax");
  
  cout<<"MKBar: "<<MKaonBar<<" GeV"<<endl;
  cout<<"MEtaExp: "<<MEtaExp<<" GeV"<<endl;
  cout<<"MEtaExp: "<<ratioExp<<endl;
  
  const djvec_t Pion=load_meson(0,0);
  const djvec_t Nucleon=load_nucleon(0);
  const djvec_t Pion_eff=effective_mass(Pion);
  const djack_t mpion=constant_fit(Pion_eff,tmin,tmax,combine("plots/pion.xmg"));
  cout<<"aMPi: "<<mpion.ave_err()<<endl;
  cout<<"aM2Pi: "<<sqr(mpion).ave_err()<<endl;
  const djvec_t Nucleon_eff=effective_mass(Nucleon);
  const djack_t mnucleon=constant_fit(Nucleon_eff,tmin,tmax,combine("plots/nucleon.xmg"));
  cout<<"aMNuc: "<<mnucleon.ave_err()<<endl;
  
  djvec_t ratio(nq),mkaon(nq),momega(nq),meta_ss(nq),ratioRec(nq);
  for(size_t iq=0;iq<nq;iq++)
    {
      const djvec_t Omega=load_omega(iq);
      const djvec_t Kaon=load_meson(0,iq);
      
      const djvec_t Omega_eff=effective_mass(Omega);
      const djvec_t Kaon_eff=effective_mass(Kaon);
      const djvec_t EtaSS_eff=sqrt(2*sqr(Kaon_eff)-sqr(Pion_eff));
      const djvec_t ratio_corr=Omega_eff/EtaSS_eff;
      
      ratio[iq]=constant_fit(ratio_corr,tmin,tmax,combine("plots/ratio_%zu.xmg",iq));
      momega[iq]=constant_fit(Omega_eff,tmin,tmax,combine("plots/omega_%zu.xmg",iq));
      mkaon[iq]=constant_fit(Kaon_eff,tmin,tmax,combine("plots/kaon_%zu.xmg",iq));
      meta_ss[iq]=constant_fit(EtaSS_eff,tmin,TH-3,combine("plots/etaSS_%zu.xmg",iq));
      ratioRec[iq]=momega[iq]/sqrt(2*sqr(mkaon[iq])-sqr(mpion));
      
      cout<<"MK: "<<mkaon[iq].ave_err()<<endl;
      cout<<"MEtaSS: "<<meta_ss[iq].ave_err()<<endl;
      cout<<"MOmega: "<<momega[iq].ave_err()<<endl;
    }
  
  /// Fit omega mass with a second order polynomial
  grace_file_t momega_plot("plots/momega.xmg");
  momega_plot.write_vec_ave_err(amq,momega.ave_err());
  const djvec_t momegaPars=poly_fit(amq,momega,2,amq[0]*1.1);
  
  grace_file_t ratio_plot("plots/ratio.xmg");
  //ratio_plot.write_vec_ave_err(amq,ratioRec.ave_err());
  for(size_t iq=1;iq<nq;iq++)
    ratio_plot.write_ave_err(amq[iq],ratio[iq].ave_err());
  
  /// Fit the ratio with a second order polynomial
  const djvec_t ratioPars=poly_fit(amq,ratio,2,amq[1]*0.9);
  ratio_plot.write_polygon([ratioPars](const double x){return poly_eval(ratioPars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::GREEN4);
  //cout<<ratioPars.ave_err()<<endl;
  const djack_t ams_phys=parab_solve(ratioPars,ratioExp);
  cout<<"ams(quad): "<<ams_phys.ave_err()<<endl;
  
  /// Fit the ratio at all ms with a first order polynomial
  const djvec_t ratioLin3Pars=poly_fit(amq,ratio,1,amq[1]*0.9);
  ratio_plot.write_polygon([ratioLin3Pars](const double x){return poly_eval(ratioLin3Pars,x);},amq[1]*0.9,amq[nq-1]*1.1,grace::VIOLET);
  //cout<<ratioLinPars.ave_err()<<endl;
  const djack_t ams_lin3_phys=-(ratioLin3Pars[0]-ratioExp)/ratioLin3Pars[1];
  cout<<"ams(lin3): "<<ams_lin3_phys.ave_err()<<endl;
  
  /// Fit the ratio at all ms with a first order polynomial
  const djvec_t ratioLin2Pars=poly_fit(amq,ratio,1,amq[2]*0.9);
  ratio_plot.write_polygon([ratioLin2Pars](const double x){return poly_eval(ratioLin2Pars,x);},amq[2]*0.9,amq[nq-1]*1.1,grace::VIOLET);
  //cout<<ratioLinPars.ave_err()<<endl;
  const djack_t ams_lin2_phys=-(ratioLin2Pars[0]-ratioExp)/ratioLin2Pars[1];
  cout<<"ams(lin2): "<<ams_lin2_phys.ave_err()<<endl;
  
  /// Finds the physical ams
  ratio_plot.set_all_colors(grace::ORANGE);
  ratio_plot.write_ave_err(ams_phys.ave_err(),{ratioExp,0.0});
  
  const djack_t momega_inte=poly_eval(momegaPars,ams_phys);
  
  const djack_t lata_inv=MOmega/momega_inte;
  const djack_t lata=1/lata_inv;
  cout<<"latinv "<<lata_inv.ave_err()<<" GeV"<<endl;
  cout<<"lat "<<lata.ave_err()<<" GeV^-1"<<endl;
  
  const double rat2_phys=MPionZero/MOmega;
  const djack_t rat2=mpion/momega_inte/rat2_phys;
  cout<<"rat2: "<<sqr(rat2).ave_err()<<endl;
  
  return 0;
}
