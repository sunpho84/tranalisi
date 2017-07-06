#include <tranalisi.hpp>

#include <fpiem_FSE.hpp>

const double a=0.457;
//const double a2=a*a;
const double MPi_phys=0.1349766,fPi_phys=0.13041;
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

//return xi
template <class T> T xi_fun(const T &mpi,const T &fpi)
{return sqr((T)(mpi/(4.0*M_PI*fpi)));}

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

const double xi_phys=xi_fun(MPi_phys,fPi_phys);

double fpi_inf_inv_fun(double si,double xii,double p6,double p1,double p2)
{
  double ell_6=p6-log(xii/xi_phys);
  double Rsi=2.0/3.0+(1.0+4.0/si)*(2.0+sqrt(1.0+4.0/si)*log((sqrt(si+4)-sqrt(si))/(sqrt(si+4)+sqrt(si))));
  double ans=1.0+si*xii*(ell_6-1.0+Rsi)/3.0+sqr(xii)*si*(p1+p2*si)/6.0;
  //cout<<"Rs: "<<Rsi<<" xi: "<<xii<<", xi_phys: "<<xi_phys<<", ans: "<<ans<<", si: "<<si<<", p6: "<<p6<<" p1: "<<p1<<", p2: "<<p2<<endl;
  return ans;
}

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
  
  vector<djvec_t> ff(nens_used),ff_FSE(nens_used);
  valarray<valarray<double>> a2Q2(nens_used);
  djvec_t aMPi(nens_used),afPi(nens_used);
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
      
      vector<double> a2p2(nth),pi(nth);
      djvec_t Z2Pi(nth),aEfit(nth),aEcont(nth),aElat(nth),aE(nth),a2P2(nth),aP0(nth);
      vector<djvec_t> corr_PP(nth);
      
      ff[iens].resize(nth);
      ff_FSE[iens].resize(nth);
      a2Q2[iens].resize(nth);
      
      for(size_t ith=0;ith<nth;ith++)
	{
	  corr_PP[ith]=load_corr("pp",ith,1,ens);
	  size_t TH=ens.T/2;
	  
	  two_pts_fit(Z2Pi[ith],aEfit[ith],corr_PP[ith],TH,ens.tmin_2pts,ens.tmax_2pts,ppath+"/pp_effmass_"+to_string(ith)+".xmg");
	  
	  //compute kinematics
	  pi[ith]=ens.th[ith]*M_PI/ens.L;
	  a2p2[ith]=3.0*sqr(pi[ith]);
	  aEcont[ith]=cont_en(aEfit[0],pi[ith]);
	  aElat[ith]=latt_en(aEfit[0],pi[ith]);
	  aE[ith]=aElat[ith];
	  
	  a2Q2[iens][ith]=4.0*a2p2[ith];
	  aP0[ith]=2.0*aE[ith];
	  a2P2[ith]=sqr(aP0[ith]);
	  
	  eff_all.write_vec_ave_err(effective_mass(corr_PP[ith]).ave_err());
	  eff_all.write_constant_band(ens.tmin_2pts,ens.tmax_2pts,aEfit[ith]);
	}
      
      //decay constant
      afPi[iens]=2*ens.aml*sqrt(Z2Pi[0])/sqr(aEfit[0]);
      aMPi[iens]=aEfit[0];
      
      //write the dispersion relation
      grace_file_t disprel(ppath+"disprel.xmg");
      disprel.write_vec_ave_err(a2p2,djvec_t(sqr(aEfit)).ave_err());
      disprel.write_polygon([&aEfit](double a2p2){return djack_t(a2p2+sqr(aEfit[0]));},0,a2p2[nth-1]);
      disprel.write_polygon([&aEfit](double a2p2){return djack_t(sqr(latt_en(aEfit[0],sqrt(a2p2/3))));},0,a2p2[nth-1]);
      
      auto aperiodic_effective_mass=[](const djvec_t corr){return forward_derivative(djvec_t(log(corr)));};
      
      djvec_t mel(nth);
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
	  ff[iens][ith]=mel[ith]*2*aE[ith]*exp(aE[ith]*tsep); //see eq.20 of 0812.4042
	  ff_FSE[iens][ith]=FSE_V(aE[0],ens.L,afPi[iens],ens.th[ith]/2);
	  //ff[ith]=mel[ith]/corr_PP[ith][tsep]; //does not work better
	}
      
      //normalize and renormalize
      ff[iens]/=(djack_t)(ff[iens][0]);
      
      //plot ff
      grace_file_t ff_plot(ppath+"ff.xmg");
      ff_plot.write_vec_ave_err(a2Q2[iens],ff[iens].ave_err());
      
      ff_all.write_vec_ave_err(a2Q2[iens],ff[iens].ave_err());
      ff_all.set_legend(ens.path);
      
      for(size_t ith=0;ith<nth;ith++)
	table<<ith<<"\t"<<a2p2[ith]<<"\t"<<a2Q2[iens][ith]<<"\t"<<aEfit[ith]<<"\t"<<ff[iens][ith]<<"\t"<<ff_FSE[iens][ith]<<endl;
      
      table<<endl;
      table<<"\tafPi="<<afPi[iens].ave_err()<<"\tML: "<<djack_t(aMPi[iens]*ens.L).ave_err()<<endl;
    }
  
  //fitting
  jack_fit_t fitter;
  djack_t C,LEC_6,B1,B2;
  djvec_t xi=xi_fun(aMPi,afPi);
  vector<djvec_t> s(nens_used);
  for(size_t iens=0;iens<nens_used;iens++)
    {
      s[iens].resize(ens_data[iens].nth());
      for(size_t ith=0;ith<ens_data[iens].nth();ith++)
	s[iens][ith]=a2Q2[iens][ith]/sqr(aMPi[iens]);
    }
  
  size_t iC=fitter.add_fit_par(C,"C",{11.9,0.1});
  size_t iB1=fitter.add_fit_par(B1,"B1",{54.3,0.1});
  size_t iB2=fitter.add_fit_par(B2,"B2",{17.9,0.1});
  size_t iLEC_6=fitter.add_fit_par(LEC_6,"LEC_6",{15.9,0.1});
  fitter.fix_par_to(iB1,0.0);
  fitter.fix_par_to(iB2,0.0);
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      size_t nth=ens.nth();
      
      for(size_t ith=1;ith<nth;ith++)
	fitter.add_point(1/ff[iens][ith],
			 [iC,iB1,iB2,iLEC_6,&xi,&s,ith,iens,&ff_FSE]
			 (const vector<double> &p,int iel)
			 {
			   double si=s[iens][ith][iel];
			   double xii=xi[iens][iel];
			   double fpi_inf_inv_ans=fpi_inf_inv_fun(si,xii,p[iLEC_6],p[iB1],p[iB2]);
			   //cout<<fpi_inf_inv_ans<<" "<<-p[iC]*ff_FSE[iens][ith][iel]*sqr(fpi_inf_inv_ans)<<endl;
			   return fpi_inf_inv_ans*(1-p[iC]*ff_FSE[iens][ith][iel]*fpi_inf_inv_ans);
			 });
    }
  
  fit_debug=1;
  fitter.fit();
  
  return 0;
}
