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
  bool use,use_for_FSE,use_for_chir;
  vector<double> th;
  size_t nth(){return th.size();}
  size_t T,L,spat_vol;
  size_t tmin_2pts,tmax_2pts;
  size_t tmin_3pts,tmax_3pts;
  double aml;
  string path;
};

//! fit range variations
vector<int> fit_range_variations={0,2};

//! NLO or NNLO fit
enum chirord_t{NLO,NNLO};
vector<chirord_t> chirord_variations={NLO,NNLO};
vector<string> chirord_tag={"NLO","NNLO"};

//! pion fit
enum pion_masses_t{PHYS_PI,ALL_PI};
vector<pion_masses_t> pion_masses_variations={PHYS_PI,ALL_PI};
vector<string> pion_masses_tag={"PHYS","ALL"};

//! no FSE, Colangelo or expanded version
enum FSE_t{NO_FSE,COLANGELO,EXPANDED};
vector<FSE_t> FSE_variations={NO_FSE,COLANGELO,EXPANDED};
vector<string> FSE_tag={"None","Colangelo","Expanded"};

//! maximal value of s
vector<double> s_max_variations={2,0};

//! use or not correlation and corrections
bool use_cov;
double cov_fact;

//! analytic or numerical way of removing interpolating fields
enum RAT_t{AN_RAT,NU_RAT};
vector<RAT_t> RAT_variations={AN_RAT,NU_RAT};
vector<string> RAT_tag={"Analytic","Numerical"};

vector<ens_data_t> ens_data;

enum{iRAT_comp,ifit_range_comp,iuse_pion_masses_comp,iuse_chirord_comp,iuse_FSE_comp,is_max_comp};
index_t syst_ind({
    {"Ratio",RAT_variations.size()},
    {"FitRange",fit_range_variations.size()},
    {"Pion masses",pion_masses_variations.size()},
    {"Chir Ord",chirord_variations.size()},
    {"FSE",FSE_variations.size()},
    {"Smax",s_max_variations.size()}});

//! return xi
template <class T>
T xi_fun(const T &mpi,const T &fpi)
{return sqr((T)(mpi/(4.0*M_PI*fpi)));}

//! return s
template <class T>
T si_fun_Q2(const double &Q2,const T &mpi)
{return Q2/sqr(mpi);}

//! return Q2
template <class T>
T Q2_fun_si(const double &si,const T &mpi)
{return si*sqr(mpi);}

//! return theta given a q2
template <class T>
T th_fun_Q2(const T &Q2,int L)
{return L*sqrt(Q2/12.0)/M_PI;}

//! return q2 given theta
template <class T>
T Q2_fun_th(const T &th,int L)
{return 12.0*sqr(th*M_PI/L);}

//! get an element out of a possible list, through a key
template <class T>
T find_key(const string &key,const map<string,T> &poss)
{
  auto it=poss.find(key);
  if(it==poss.end())
    {
      cerr<<"List of possible keys: "<<endl;
      for(auto &p : poss) cerr<<p.first<<endl;
      CRASH("Unable to fine key: %s",key.c_str());
    }
  return it->second;
}

//! initialize the program
inline void fpiem_initialize(int narg,char **arg)
{
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  size_t ext_njacks=input.read<size_t>("NJacks");
  set_njacks(ext_njacks);
  
  cout.precision(16);
  size_t nens=input.read<size_t>("NEnsemble");
  
  use_cov=input.read<bool>("UseCov");
  cov_fact=input.read<double>("CovFact");
  
  input.expect({"Use","UseChir","UseFSE","L","T","t2pts","t3pts","aml","path","nth"});
  ens_data.resize(nens);
  for(size_t iens=0;iens<nens;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      
      input.read(ens.use);
      input.read(ens.use_for_chir);
      input.read(ens.use_for_FSE);
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

template <class T>
T FSE_fun(const T &mpi,double L,const T&fpi,double thhalf,FSE_t use_FSE)
{
  T out;
  switch(use_FSE)
  {
  case NO_FSE: out=0.0;break;
  case COLANGELO: out=FSE_V(mpi,L,fpi,thhalf);break;
  case EXPANDED:
    T xi=xi_fun(mpi,fpi);
    out=xi*exp(-mpi*L)/(mpi*L,1.5);
  }
  
  return out;
}

//! fit ansatz for the inverse
template <class Tpars,class Tx>
Tpars fpi_inf_inv_fun(double a2Q2,Tx aMPi,Tx afPi,Tx ff_FSE,Tpars p6,Tpars p1,Tpars p2,Tpars pC1,Tpars pC2,FSE_t use_FSE)
{
  Tx si=si_fun_Q2(a2Q2,aMPi);
  Tx xi=xi_fun(aMPi,afPi);
  Tpars ell_6=p6-log(xi/xi_phys);
  Tpars Rsi=2.0/3.0+(1.0+4.0/si)*(2.0+sqrt(1.0+4.0/si)*log((sqrt(si+4)-sqrt(si))/(sqrt(si+4)+sqrt(si))));
  Tpars ans=1.0+si*xi*(ell_6-1.0+Rsi)/3.0+sqr(xi)*si*(p1+p2*si)/6.0;
  Tpars FSE_prefact=0.0;
  switch(use_FSE)
    {
    case NO_FSE: FSE_prefact=0;break;
    case COLANGELO: FSE_prefact=pC1;break;
    case EXPANDED:  FSE_prefact=pC1*si+pC2*si*si;break;
    }
  Tpars FSE_fact=(1-FSE_prefact*ff_FSE*ans);
  //cout<<"Rs: "<<Rsi<<" xi: "<<xii<<", xi_phys: "<<xi_phys<<", ans: "<<ans<<", si: "<<si<<", p6: "<<p6<<" p1: "<<p1<<", p2: "<<p2<<endl;
  return ans*FSE_fact;
}

//! fitting
djack_t fit_fpiinv(const djvec_t &aMPi,const djvec_t &afPi,const valarray<valarray<double>> &a2Q2,const vector<djvec_t> &ff,const vector<djvec_t> &ff_FSE,pion_masses_t use_pion_masses,chirord_t use_chirord,FSE_t use_FSE,double s_max,size_t isyst,bool cov_flag=false,double eps_cov_flag=0.0)
{
  jack_fit_t fitter;
  djack_t C1,C2,LEC_6,B1,B2;
  size_t iC1=fitter.add_fit_par(C1,"C1",{11.9*0,0.1});
  size_t iC2=fitter.add_fit_par(C2,"C2",{11.9*0,0.1});
  size_t iB1=fitter.add_fit_par(B1,"B1",{54.3,0.1});
  size_t iB2=fitter.add_fit_par(B2,"B2",{17.9,0.1});
  size_t iLEC_6=fitter.add_fit_par(LEC_6,"LEC_6",{15.9,0.1});
  
  switch(use_pion_masses)
    {
    case PHYS_PI:
      fitter.fix_par_to(iB1,0.0);
      break;
    case ALL_PI:
      break;
    }
  
  switch(use_chirord)
    {
    case NLO:
      fitter.fix_par_to(iB1,0.0);
      fitter.fix_par_to(iB2,0.0);
      break;
    case NNLO:
      break;
    }
  
  switch(use_FSE)
    {
    case NO_FSE:
      fitter.fix_par_to(iC1,0.0);
      fitter.fix_par_to(iC2,0.0);
      break;
    case COLANGELO:
      fitter.fix_par_to(iC2,0.0);
      break;
    case EXPANDED:
      break;
    }
  
  vector<bool> used(ens_data.size());
  for(size_t iens=0;iens<ens_data.size();iens++)
    {
      ens_data_t &ens=ens_data[iens];
      size_t nth=ens.nth();
      used[iens]=(ens.use and (use_FSE!=NO_FSE or ens.use_for_FSE) and (use_pion_masses==ALL_PI or ens.use_for_chir));
      
      if(used[iens])
	for(size_t ith=1;ith<nth;ith++)
	  if((s_max==0 or si_fun_Q2(a2Q2[iens][ith],aMPi[iens].ave())<s_max))
	    fitter.add_point(1/ff[iens][ith],
			     [&a2Q2,iC1,iC2,iB1,iB2,iLEC_6,&aMPi,&afPi,ith,iens,&ff_FSE,&use_FSE]
			     (const vector<double> &p,int iel)
			     {return fpi_inf_inv_fun(a2Q2[iens][ith],aMPi[iens][iel],afPi[iens][iel],ff_FSE[iens][ith][iel],p[iLEC_6],p[iB1],p[iB2],p[iC1],p[iC2],use_FSE);}
			     ,iens);
    }
  
  //fit_debug=1;
  fitter.fit(cov_flag,eps_cov_flag);
  
  //write plots
  for(size_t iens=0;iens<ens_data.size();iens++)
    {
      vector<size_t> syst_comp=syst_ind(isyst);
      size_t iRAT=syst_comp[iRAT_comp];
      size_t ifit_range=syst_comp[ifit_range_comp];
      size_t iuse_pion_masses=syst_comp[iuse_pion_masses_comp];
      size_t iuse_chirord=syst_comp[iuse_chirord_comp];
      size_t iuse_FSE=syst_comp[iuse_FSE_comp];
      int dtf=fit_range_variations[ifit_range];
      
      grace_file_t plot_mfix("plots/fit"+to_string(isyst)+"_inv_fpi_fun_q2_ens"+to_string(iens)+".xmg");
      plot_mfix.set_title(ens_data[iens].path+", "+(used[iens]?"":"not ")+"used");
      plot_mfix.set_subtitle_size(1.0);
      plot_mfix.set_subtitle("Ratio: "+RAT_tag[iRAT]+", "
			     "FSE: "+FSE_tag[iuse_FSE]+", "
			     "pion masses: "+pion_masses_tag[iuse_pion_masses]+", "
			     "chirord: "+chirord_tag[iuse_chirord]+", "
			     "dt fit: "+to_string(dtf)+", "
			     "s max: "+to_string(s_max));
      plot_mfix.set_xaxis_label("$$a^2Q^2");
      plot_mfix.set_yaxis_label("$$1/f_+^\\pi");
      
      //write data, in 2 scansions
      djvec_t inv=1/ff[iens];
      plot_mfix.set_settype(grace::XYDY);
      for(size_t it=0;it<2;it++) //0=used, 1=not used
	{
	  if(it) plot_mfix.new_data_set();
	  
	  if(it==1) plot_mfix.set_symbol_fill_pattern(grace::FILLED_SYMBOL);
	  for(size_t ith=1;ith<ens_data[iens].nth();ith++)
	    {
	      bool th_used=(used[iens] and (s_max==0 or si_fun_Q2(a2Q2[iens][ith],aMPi[iens].ave())<s_max));
	      if(it!=th_used) plot_mfix<<a2Q2[iens][ith]<<" " <<inv[ith].ave_err()<<endl;
	    }
	  plot_mfix.new_data_set();
	}
      
      size_t L=ens_data[iens].L;
      //select the range of Q2 to paint in the two colors
      double Q2_min=a2Q2[iens][1]/2;
      double Q2_max=a2Q2[iens][ens_data[iens].nth()-1];
      double Q2_int=Q2_fun_si(s_max,aMPi[iens]).ave();
      if(Q2_int<Q2_min) Q2_int=Q2_min;
      if(Q2_int>Q2_max) Q2_int=Q2_max;
      double xmin[2]={Q2_min,Q2_int};
      double xmax[2]={Q2_int,Q2_max};
      for(size_t it=0;it<2;it++)
	plot_mfix.write_polygon([&L,&aMPi,afPi,iens,LEC_6,B1,B2,C1,C2,use_FSE,xmin,xmax,it](double x)
				{
				  auto th=th_fun_Q2(x,L);
				  auto F=FSE_fun(aMPi[iens],L,afPi[iens],th/2.0,use_FSE);
				  //cout<<th<<" "<<F<<endl;
				  return fpi_inf_inv_fun(x,aMPi[iens],afPi[iens],F,LEC_6,B1,B2,C1,C2,use_FSE);}
				,xmin[it],xmax[it],50);
    }
  
  cout<<"C1: "<<C1<<endl;
  cout<<"C2: "<<C2<<endl;
  cout<<"LEC_6: "<<LEC_6<<endl;
  cout<<"B1: "<<B1<<endl;
  cout<<"B2: "<<B2<<endl;
  
  double fpi_exp=0.1304;
  djack_t chrad=1.0/sqr(4*M_PI*fpi_exp)*(2*(LEC_6-1.0)+B1*xi_phys)*sqr(0.19731);
  cout<<"Charge radius: "<<chrad<<" fm^2"<<endl;
  
  return chrad;
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
  
  cout<<"Studying "<<syst_ind.max()<<" systematics"<<endl;
  
  djvec_t ch_rad(syst_ind.max());
  for(size_t isyst=0;isyst<syst_ind.max();isyst++)
    {
      cout<<"################################################## "<<isyst<<" ##################################################"<<endl;
      
      vector<size_t> syst_comp=syst_ind(isyst);
      size_t iRAT=syst_comp[iRAT_comp];
      size_t ifit_range=syst_comp[ifit_range_comp];
      size_t iuse_pion_masses=syst_comp[iuse_pion_masses_comp];
      size_t iuse_chirord=syst_comp[iuse_chirord_comp];
      size_t iuse_FSE=syst_comp[iuse_FSE_comp];
      size_t is_max=syst_comp[is_max_comp];
      
      FSE_t use_FSE=FSE_variations[iuse_FSE];
      pion_masses_t use_pion_masses=pion_masses_variations[iuse_pion_masses];
      chirord_t use_chirord=chirord_variations[iuse_chirord];
      int dtf=fit_range_variations[ifit_range];
      double s_max=s_max_variations[is_max];
      
      cout<<"Ratio: "<<RAT_tag[iRAT]<<endl;
      cout<<"Use FSE: "<<FSE_tag[iuse_FSE]<<endl;
      cout<<"Use pion masses: "<<pion_masses_tag[iuse_pion_masses]<<endl;
      cout<<"Use chirord: "<<chirord_tag[iuse_chirord]<<endl;
      cout<<"dt for fit: "<<dtf<<endl;
      cout<<"s max: "<<s_max<<endl;
      
      vector<djvec_t> ff(ens_data.size()),ff_FSE(ens_data.size());
      valarray<valarray<double>> a2Q2(ens_data.size());
      djvec_t aMPi(ens_data.size()),afPi(ens_data.size());
      for(size_t iens=0;iens<ens_data.size();iens++)
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
	      
	      two_pts_fit(Z2Pi[ith],aEfit[ith],corr_PP[ith],TH,ens.tmin_2pts+dtf,ens.tmax_2pts-dtf,ppath+"/pp_effmass_"+to_string(ith)+".xmg");
	      
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
	      eff_all.write_constant_band(ens.tmin_2pts+dtf,ens.tmax_2pts-dtf,aEfit[ith]);
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
	      mel[ith]=constant_fit(vector_ff,ens.tmin_3pts+dtf,ens.tmax_3pts-dtf,ens.path+"/plots/vv_th"+to_string(ith)+".xmg");
	      djvec_t vector_ff_effmass=aperiodic_effective_mass(vector_ff);
	      vector_ff_effmass.ave_err().write(ens.path+"/plots/vv_effmass_th"+to_string(ith)+".xmg");
	      
	      vv_all.write_vec_ave_err(vector_ff.ave_err());
	      vv_all.write_constant_band(ens.tmin_3pts+dtf,ens.tmax_3pts-dtf,mel[ith]);
	      
	      //compute ff
	      size_t tsep=ens.T/2;
	      switch(iRAT)
		{
		case AN_RAT: ff[iens][ith]=mel[ith]*2*aE[ith]*exp(aE[ith]*tsep);break; //see eq.20 of 0812.4042
		case NU_RAT: ff[iens][ith]=mel[ith]/corr_PP[ith][tsep];break; //does not work better
		}
	      
	      ff_FSE[iens][ith]=FSE_fun(aMPi[iens],ens.L,afPi[iens],ens.th[ith]/2,use_FSE);
	      //cout<<"check:"<<ens.th[ith]<<" "<<ff_FSE[iens][ith]<<endl;
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
      
      ch_rad[isyst]=fit_fpiinv(aMPi,afPi,a2Q2,ff,ff_FSE,use_pion_masses,use_chirord,use_FSE,s_max,isyst,use_cov,cov_fact);
    }
  
  perform_analysis(ch_rad,syst_ind,"Charge radius");
  
  return 0;
}
