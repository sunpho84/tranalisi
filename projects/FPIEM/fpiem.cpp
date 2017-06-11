#include <tranalisi.hpp>

#include <fpiem_FSE.hpp>

const double a=0.457,a2=a*a;
size_t EVN=1,ODD=-1;

//! hold the data for a single ensemble
class ens_data_t
{
public:
  vector<double> th;
  size_t nth(){return th.size();}
  size_t T,L,spat_vol;
  size_t tmin_2pts,tmax_2pts;
  size_t tmin_3pts,tmax_3pts;
  double aml;
  string path;
};

vector<ens_data_t> ens_data;
size_t nens_used;

//! initialize the program
inline void fpiem_initialize(int narg,char **arg)
{
  set_njacks(15);
  
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  nens_used=input.read<int>("NEnsemble");
  
  input.expect({"L","T","t2pts","t3pts","aml","path","nth"});
  ens_data.resize(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      
      input.read(ens.L);
      ens.spat_vol=ens.L*ens.L*ens.L;
      input.read(ens.T);
      input.read(ens.tmin_2pts);
      input.read(ens.tmax_2pts);
      input.read(ens.tmin_3pts);
      input.read(ens.tmax_3pts);
      input.read(ens.aml);
      input.read(ens.path);
      size_t nth;
      input.read(nth);
      ens.th.resize(nth);
      for(size_t ith=0;ith<nth;ith++) input.read(ens.th[ith]);
    }
}

//! load the 2pts
djvec_t load_corr(const char *name,size_t ith,int par,const ens_data_t &ens)
{return read_djvec(ens.path+"/jacks/"+name+"."+combine("%02zu",ith),ens.T).symmetrized(par);}

int main(int narg,char **arg)
{
  fpiem_initialize(narg,arg);
  
  grace_file_t ff_all("plots/ff_all.xmg");
  ff_all.set_xaxis_label("-(aQ)\\S2\\N/"+combine("(%lg GeV\\S-1\\N)",a)+"\\S2\\N");
  
  string table_path="tables/all.txt";
  ofstream table(table_path);
#define SEP "\t   "
  table<<"Ith" SEP "(ap)^2" SEP "(aQ)^2" SEP " aE\t" SEP "ff"<<endl;
  
  if(not table.good()) CRASH("Opening %s",table_path.c_str());
  table<<fixed;
  table.precision(8);
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      string ppath=ens.path+"/plots/";
      size_t nth=ens.nth();
      
      table<<endl;
      table<<"-----------------------------"<<ens.path<<"---------------------------"<<endl;
      
      grace_file_t eff_all(ppath+"/eff_all.xmg");
      grace_file_t vv_all(ppath+"/vv_all.xmg");
      vv_all.set_title("<\\xp\\0(p)|V\\s0\\N|\\xp\\0(-p)>");
      
      for(auto &it : {&eff_all,&vv_all})
	{
	  it->set_subtitle(ens.path);
	  auto color_scheme={grace::RED,grace::GREEN4,grace::ORANGE,grace::VIOLET,grace::BLUE,grace::BROWN,grace::MAGENTA};
	  it->set_line_color_scheme(color_scheme);
	  it->set_color_scheme(color_scheme);
	}
      
      vector<double> a2p2(nth),a2Q2(nth),Q2(nth),pi(nth);
      djvec_t Z2Pi(nth),aEfit(nth),aEcont(nth),aElat(nth),aE(nth),a2P2(nth),aP0(nth);
      djack_t afPi;
      vector<djvec_t> corr_PP(nth);
      
      for(size_t ith=0;ith<nth;ith++)
	{
	  corr_PP[ith]=load_corr("pp",ith,1,ens);
	  size_t TH=ens.T/2;
	  
	  two_pts_fit(Z2Pi[ith],aEfit[ith],corr_PP[ith],TH,ens.tmin_2pts,ens.tmax_2pts,ppath+"/pp_effmass_"+to_string(ith)+".xmg");
	  
	  //compute kinematics
	  pi[ith]=ens.th[ith]*M_PI/ens.L;
	  a2p2[ith]=3*sqr(pi[ith]);
	  aEcont[ith]=cont_en(aEfit[0],pi[ith]);
	  aElat[ith]=latt_en(aEfit[0],pi[ith]);
	  aE[ith]=aElat[ith];
	  
	  a2Q2[ith]=4*a2p2[ith];
	  Q2[ith]=a2Q2[ith]/a2;
	  aP0[ith]=2*aE[ith];
	  a2P2[ith]=sqr(aP0[ith]);
	  
	  eff_all.write_vec_ave_err(effective_mass(corr_PP[ith]).ave_err());
	  eff_all.write_constant_band(ens.tmin_2pts,ens.tmax_2pts,aEfit[ith]);
	}
      
      //decay constant
      afPi=2*ens.aml*sqrt(Z2Pi[0])/sqr(aEfit[0]);
      
      //write the dispersion relation
      grace_file_t disprel(ppath+"disprel.xmg");
      disprel.write_vec_ave_err(a2p2,djvec_t(sqr(aEfit)).ave_err());
      disprel.write_polygon([&aEfit](double a2p2){return djack_t(a2p2+sqr(aEfit[0]));},0,a2p2[nth-1]);
      disprel.write_polygon([&aEfit](double a2p2){return djack_t(sqr(latt_en(aEfit[0],sqrt(a2p2/3))));},0,a2p2[nth-1]);
      
      auto aperiodic_effective_mass=[](const djvec_t corr){return forward_derivative(djvec_t(log(corr)));};
      
      djvec_t mel(nth),ff(nth);
      for(size_t ith=0;ith<nth;ith++)
	{
	  //load and improve
	  djvec_t vector_ff=load_corr("vv",ith,ODD,ens);
	  vector_ff=djvec_t(vector_ff+vector_ff.inverse()).subset(0,ens.T/4+1)/2;
	  vector_ff[0]=vector_ff[1];
	  
	  //extract matrix element and print effective mass
	  mel[ith]=constant_fit(vector_ff,ens.tmin_3pts,ens.tmax_3pts,ens.path+"/plots/vv_th"+to_string(ith)+".xmg");
	  djvec_t vector_ff_effmass=aperiodic_effective_mass(vector_ff);
	  vector_ff_effmass.ave_err().write(ens.path+"/plots/vv_effmass_th"+to_string(ith)+".xmg");
	  
	  vv_all.write_vec_ave_err(vector_ff.ave_err());
	  vv_all.write_constant_band(ens.tmin_3pts,ens.tmax_3pts,mel[ith]);
	  
	  //compute ff
	  size_t tsep=ens.T/2;
	  ff[ith]=mel[ith]*2*aE[ith]*exp(aE[ith]*tsep); //see eq.20 of 0812.4042
	  
	  //ff[ith]=mel[ith]/corr_PP[ith][tsep]; //does not work better
	}
      
      //normalize and renormalize
      ff/=(djack_t)(ff[0]);
      
      //plot ff
      grace_file_t ff_plot(ppath+"ff.xmg");
      ff_plot.write_vec_ave_err(Q2,ff.ave_err());
      
      ff_all.write_vec_ave_err(Q2,ff.ave_err());
      ff_all.set_legend(ens.path);
      
      for(size_t ith=0;ith<nth;ith++)
	table<<ith<<"\t"<<a2p2[ith]<<"\t"<<a2Q2[ith]<<"\t"<<aEfit[ith].ave_err()<<"\t"<<ff[ith].ave_err()<<endl;
      
      table<<endl;
      table<<"\tafPi="<<afPi.ave_err()<<"\tML: "<<djack_t(aEfit[0]*ens.L).ave_err()<<endl;
      //"\tfPi="<<smart_print(djack_t(afPi/a).ave_err())<<endl;
    }
  
  return 0;
}
