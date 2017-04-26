
#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <Kl2_IB_FSE.hpp>
#include <common.hpp>

////////////////////////////////////////////////// quarks /////////////////////////////////////////////////

class qpars_t
{
public:
  size_t iq;
  double dm;
  double eq;
  qpars_t(size_t iq,double dm,double eq) : iq(iq),dm(dm),eq(eq) {}
};
const qpars_t qU(0,1,eu),qD(0,-1,ed),qS(1,0,es),qC(2,0,ec);
const vector<qpars_t> all_qpars({qU,qD,qS,qC});

////////////////////////////////////////////////// QCD mesons //////////////////////////////////////////////

//! hold name and consittuent quark for a given meson
class QCD_mes_pars_t
{
public:
  size_t iq1; //< index of quark 1
  size_t iq2; //< index of quark 2
  size_t itint; //< which interval range this meson has to use (0=Pi, 1=K, 2=D,Ds)
  string name; //< name of the quark
  QCD_mes_pars_t(const qpars_t &q1,const qpars_t &q2,size_t itint,const string &name) : iq1(q1.iq),iq2(q2.iq),itint(itint),name(name) {}
};

const vector<QCD_mes_pars_t> QCD_mes_pars({{qD,qD,1,"Pi"},{qD,qS,1,"K"},{qC,qD,2,"D"},{qC,qS,2,"Ds"}});
enum{iPi,iK,iD};
const size_t &nQCD_mes=QCD_mes_pars.size();

/////////////////////////////////////////////////// QED mesons /////////////////////////////////////

//! hold name and consittuent quark for a given QED meson
class QED_mes_pars_t
{
public:
  size_t iq1; //< index of quark 1
  size_t iq2; //< index of quark 2
  double dm1; //< contribution by which scalar insertion on quark 1 mast be added
  double dm2; //< contribution by which scalar insertion on quark 2 mast be added
  double eq1; //< charge of quark 1
  double eq2; //< charge of quark 2
  size_t irev; //< quark line to be reversed
  size_t iQCD; //< index of pure QCD meson corresponding to this meson
  string name; //< name of the meson
  QED_mes_pars_t(const qpars_t &q1,const qpars_t &q2,size_t irev,size_t iQCD,const string &name) :
    iq1(q1.iq),iq2(q2.iq),dm1(q1.dm),dm2(q2.dm),eq1(q1.eq),eq2(q2.eq),irev(irev),iQCD(iQCD),name(name) {}
};

const size_t QREV1=0; //< value to revert quark 1
const size_t QREV2=1; //< value to revert quark 2
const vector<QED_mes_pars_t> QED_mes_pars({{qU,qD,QREV2,0,"PiPlus"},
					   {qS,qU,QREV1,1,"KPlus"},
					   {qS,qD,QREV1,1,"K0"},
					   {qC,qD,QREV2,2,"DPlus"},
					   {qC,qU,QREV2,2,"D0"},
					   {qC,qS,QREV2,3,"DsPlus"}});
enum{iPiPlus,iKPlus,iK0,iDPlus,iD0,iDsPlus};
const size_t &nQED_mes=QED_mes_pars.size();

////////////////////////////////////////////////////////////////////

const size_t nleps=2; //< number of leptons
const size_t nmes_tint=3; //< number of time intervals

class ens_pars_t
{
public:
  size_t iult; //< input in the ultimate file
  size_t ib; //< beta index
  size_t T; //< time extent
  size_t L; //< spatial size
  double aml; //< bare light quark mass
  double MLep[nleps]; //< mass of leptons
  double MMes[4]; //< mass of mesons (0=Pi, 1=K, 2=D, 3=Ds)
  int use_for_L; //< use for FSE analysis
  string path; //< path (name)
  
  vector<size_t> tmin,tmax; //< range of fit
  ens_pars_t() : tmin(nmes_tint),tmax(nmes_tint) {}
};
vector<ens_pars_t> ens_pars; //< parameters of all ensemble
size_t nens_used; //< number of ensemble used

const size_t nqmass=3; //< number of quark mass
const size_t nr=2; //< number of r
const index_t ind_2pts({{"NMass",nqmass},{"NMass",nqmass},{"Nr",nr},{"RI",2}});

vector<size_t> iQED_mes_of_proc ({0,1,3,3,5,5}); //< index of the QED meson corresponding to a given process
vector<size_t> imlep_of_proc({0,0,0,1,0,1}); //< index of the lepton corresponding to a given process

//! initialize everything
void initialize(int narg,char **arg)
{
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  //read where to read input and how many ensemble
  string ultimate_path=input.read<string>("UltimatePath");
  init_common_IB(ultimate_path);
  nens_used=input.read<int>("NEnsemble");
  ens_pars.resize(nens_used);
  Za.resize(nbeta);
  
  input.expect({"Ens","beta","aml","L","UseForL","T","TPi","TK","TD","path","MLep0","MLep1","MMes0","MMes1","MMes2","MMes3"});
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      
      input.read(ens.iult);
      input.read(ens.ib);
      input.read(ens.aml);
      input.read(ens.L);
      input.read(ens.use_for_L);
      input.read(ens.T);
      for(size_t itint=0;itint<nmes_tint;itint++)
	{
	  input.read(ens.tmin[itint]);
	  input.read(ens.tmax[itint]);
	}
      input.read(ens.path);
      for(size_t ilep=0;ilep<nleps;ilep++) input.read(ens.MLep[ilep]);
      for(size_t iQCD_mes=0;iQCD_mes<nQCD_mes;iQCD_mes++) input.read(ens.MMes[iQCD_mes]);
    }
}

//! read a single vector, for a specific mass and r, real or imaginary
djvec_t read(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,size_t r,size_t reim)
{
  string path=combine("%s/data/corr%s",ens.path.c_str(),what);
  djvec_t out=read_djvec(path,ens.T,ind_2pts({iq1,iq2,r,reim}));
  return out;
}

//! read a combination of r and return appropriately symmetrized
djvec_t read(const char *what,const ens_pars_t &ens,int tpar,size_t iq1,size_t iq2,int rpar,size_t reim)
{
  djvec_t o(ens.T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,ens,iq1,iq2,r,reim)*(r==0?1:rpar);
  return o.symmetrized(tpar)/(1+abs(rpar));
}

//! read PP
djvec_t read_PP(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim)
{return read(combine("%s_%s",what,"P5P5").c_str(),ens,1,iq1,iq2,rpar,reim);}

//! read VP
djvec_t read_VP(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim)
{return read(combine("%s_%s",what,"V0P5").c_str(),ens,-1,iq1,iq2,rpar,reim);}

//! read AP
djvec_t read_AP(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim)
{return read(combine("%s_%s",what,"A0P5").c_str(),ens,-1,iq1,iq2,rpar,reim);}

//! compute the critical deltam
djack_t compute_deltam_cr(const ens_pars_t &ens,size_t iq,size_t iQCD_mes)
{
  string ens_qpath=ens.path+"/plots_mcrit";
  
  djvec_t V0P5_LL=read_VP("LL",ens,iq,iq,-1,IM);
  djvec_t V0P5_0M=read_VP("0M",ens,iq,iq,-1,IM);
  djvec_t V0P5_0T=read_VP("0T",ens,iq,iq,-1,IM);
  djvec_t num_deltam_cr=forward_derivative(djvec_t(V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write(combine("%s/num_deltam_cr.xmg",ens_qpath.c_str()));
  
  djvec_t V0P5_0P=read_VP("0P",ens,iq,iq,+1,RE);
  djvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write(combine("%s/den_deltam_cr.xmg",ens_qpath.c_str()));
  
  const size_t itint=QCD_mes_pars[iQCD_mes].itint;
  djack_t deltam_cr=constant_fit(djvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),ens.tmin[itint],ens.tmax[itint],combine("%s/deltam_cr_t.xmg",ens_qpath.c_str()));
  
  return deltam_cr;
}

//! read QED corrections
djvec_t read_QED(const ens_pars_t &ens,size_t iQED_mes,const djack_t &deltam_cr,const djvec_t &c_LO,
		 djvec_t(*read)(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim),const char *name)
{
  string ens_qpath=ens.path+"/plots_"+QED_mes_pars[iQED_mes].name;
  
  size_t iq1=QED_mes_pars[iQED_mes].iq1;
  size_t iq2=QED_mes_pars[iQED_mes].iq2;
  double eq1=QED_mes_pars[iQED_mes].eq1;
  double eq2=QED_mes_pars[iQED_mes].eq2;
  djvec_t c_0T1=read("0T",ens,iq2,iq1,1,RE)*eq1*eq1;
  djvec_t c_0T2=read("0T",ens,iq1,iq2,1,RE)*eq2*eq2;
  djvec_t c_0M1=read("0M",ens,iq2,iq1,1,RE)*eq1*eq1;
  djvec_t c_0M2=read("0M",ens,iq1,iq2,1,RE)*eq2*eq2;
  djvec_t c_LL=read("LL",ens,iq1,iq2,1,RE)*eq1*eq2;
  djvec_t(c_0T1/c_LO).ave_err().write(combine("%s/%s_0T1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0T2/c_LO).ave_err().write(combine("%s/%s_0T2.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0M1/c_LO).ave_err().write(combine("%s/%s_0M1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0M2/c_LO).ave_err().write(combine("%s/%s_0M2.xmg",ens_qpath.c_str(),name));
  djvec_t(c_LL/c_LO).ave_err().write(combine("%s/%s_LL.xmg",ens_qpath.c_str(),name));
  djvec_t c=-(c_LL+c_0T1+c_0M1+c_0T2+c_0M2); //minus due to slope definition
  
  djvec_t c_0P1=-(read("0P",ens,iq2,iq1,-1,IM)*eq1*eq1); //remind that i missing means
  djvec_t c_0P2=-(read("0P",ens,iq1,iq2,-1,IM)*eq2*eq2); // to take imag part changed of sign
  djvec_t(c_0P1/c_LO).ave_err().write(combine("%s/%s_0P1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0P2/c_LO).ave_err().write(combine("%s/%s_0P2.xmg",ens_qpath.c_str(),name));
  djvec_t d=djack_t(-deltam_cr)*(c_0P1+c_0P2); //minus coming from slopes
  
  return c+d;
}

//! read MASS corrections
djvec_t read_MASS(const ens_pars_t &ens,size_t iQED_mes,const djvec_t &c_LO,
		  djvec_t(*read)(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim),const char *name)
{
  string ens_qpath=ens.path+"/plots_"+QED_mes_pars[iQED_mes].name;
  const QED_mes_pars_t &pars=QED_mes_pars[iQED_mes];
  size_t iq1=pars.iq1;
  size_t iq2=pars.iq2;
  djvec_t c_0S1=read("0S",ens,iq2,iq1,1,RE)*pars.dm1;
  djvec_t c_0S2=read("0S",ens,iq1,iq2,1,RE)*pars.dm2;
  
  djvec_t(c_0S1/c_LO).ave_err().write(combine("%s/%s_0S1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0S2/c_LO).ave_err().write(combine("%s/%s_0S2.xmg",ens_qpath.c_str(),name));
  
  return -(c_0S1+c_0S2);
}

//slopes and Z for all mesons and masses
djvec_t ZP,ZA;
djvec_t M;
djvec_t DZA_QED_rel,DZA_MASS_rel;
djvec_t SL_PP_QED,SL_PP_MASS;
djvec_t SL_AP_QED,SL_AP_MASS;

index_t ind_ens_QCD_mes;
index_t ind_ens_QED_mes;
djvec_t deltam_cr;
void compute_basic_slopes()
{
  ind_ens_QCD_mes.set_ranges({{"Ens",nens_used},{"Mes",nQCD_mes}});
  ind_ens_QED_mes.set_ranges({{"Ens",nens_used},{"Mes",nQED_mes}});
  size_t nens_QCD_mes=ind_ens_QCD_mes.max();
  size_t nens_QED_mes=ind_ens_QED_mes.max();
  
  deltam_cr.resize(nens_used);
  ZP.resize(nens_QCD_mes);
  ZA.resize(nens_QCD_mes);
  M.resize(nens_QCD_mes);
  //
  DZA_QED_rel.resize(nens_QED_mes);
  DZA_MASS_rel.resize(nens_QED_mes);
  //
  SL_PP_QED.resize(nens_QED_mes);
  SL_PP_MASS.resize(nens_QED_mes);
  //
  SL_AP_QED.resize(nens_QED_mes);
  SL_AP_MASS.resize(nens_QED_mes);
  
  vector<djvec_t> jPP_LO(nens_QCD_mes);
  vector<djvec_t> jPP_MASS(nens_QED_mes);
  vector<djvec_t> jPP_QED(nens_QED_mes);
  //
  vector<djvec_t> jAP_LO(nens_QCD_mes);
  vector<djvec_t> jAP_QED(nens_QED_mes);
  vector<djvec_t> jAP_MASS(nens_QED_mes);
  
  //load everything
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      deltam_cr[iens]=compute_deltam_cr(ens,ilight,iPi);
      size_t TH=ens.T/2;
      
      //load LO for PP
      for(size_t iQED_mes=0;iQED_mes<nQED_mes;iQED_mes++)
	{
	  size_t iQCD_mes=QED_mes_pars[iQED_mes].iQCD;
	  const size_t ind_QCD=ind_ens_QCD_mes({iens,iQCD_mes});
	  const size_t ind_QED=ind_ens_QED_mes({iens,iQED_mes});
	  const size_t itint=QCD_mes_pars[iQCD_mes].itint;
	  size_t tmin=ens.tmin[itint];
	  size_t tmax=min(TH-1,ens.tmax[itint]);
	  
	  string plots_path=ens.path+"/plots_"+QED_mes_pars[iQED_mes].name;
	  
	  size_t iq1=QED_mes_pars[iQED_mes].iq1;
	  size_t iq2=QED_mes_pars[iQED_mes].iq2;
	  //
	  jPP_LO[ind_QCD]=read_PP("00",ens,iq1,iq2,1,RE);
	  jPP_MASS[ind_QED]=read_MASS(ens,iQED_mes,jPP_LO[ind_QCD],read_PP,"PP");
	  jPP_QED[ind_QED]=read_QED(ens,iQED_mes,deltam_cr[iens],jPP_LO[ind_QCD],read_PP,"PP");
	  //
	  jAP_LO[ind_QCD]=read_AP("00",ens,iq1,iq2,1,RE);
	  jAP_LO[ind_QCD].ave_err().write(plots_path+"/corr_AP_LO.xmg");
	  jAP_MASS[ind_QED]=read_MASS(ens,iQED_mes,jAP_LO[ind_QCD],read_AP,"AP");
	  jAP_QED[ind_QED]=read_QED(ens,iQED_mes,deltam_cr[iens],jAP_LO[ind_QCD],read_AP,"AP");
	  
	  djack_t ZAP,ZPP;
	  djack_t DZ_AP_QED_rel,DZ_PP_QED_rel;
	  djack_t DZ_AP_MASS_rel,DZ_PP_MASS_rel;
	  two_pts_with_ins_ratio_fit(ZAP         ,M[ind_QCD],DZ_AP_QED_rel,SL_AP_QED[ind_QED],jAP_LO[ind_QCD],jAP_QED[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_AP_LO1.xmg",plots_path+"/slope_AP_QED.xmg",-1);
	  two_pts_with_ins_ratio_fit(ZAP         ,M[ind_QCD],DZ_AP_MASS_rel,SL_AP_MASS[ind_QED],jAP_LO[ind_QCD],jAP_MASS[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_AP_LO2.xmg",plots_path+"/slope_AP_MASS.xmg",-1);
	  //
	  two_pts_with_ins_ratio_fit(ZPP,M[ind_QCD],DZ_PP_QED_rel,SL_PP_QED[ind_QED],jPP_LO[ind_QCD],jPP_QED[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_PP_LO1.xmg",plots_path+"/slope_PP_QED.xmg");
	  two_pts_with_ins_ratio_fit(ZPP,M[ind_QCD],DZ_PP_MASS_rel,SL_PP_MASS[ind_QED],jPP_LO[ind_QCD],jPP_MASS[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_PP_LO2.xmg",plots_path+"/slope_PP_MASS.xmg");
	  //
	  ZP[ind_QCD]=sqrt(ZPP);
	  ZA[ind_QCD]=ZAP/ZP[ind_QCD];
	  //
	  DZA_MASS_rel[ind_QED]=DZ_AP_MASS_rel-DZ_PP_MASS_rel/2;
	  DZA_QED_rel[ind_QED]=DZ_AP_QED_rel-DZ_PP_QED_rel/2;
	  cout<<plots_path<<", (Z)A: "<<ZA[ind_QCD].ave_err()<<endl;
	  cout<<plots_path<<", M_AP: "<<M[ind_QCD].ave_err()<<endl;
	  cout<<plots_path<<", SL_AP: "<<djack_t(SL_PP_QED[ind_QED]).ave_err()<<endl;
	  cout<<plots_path<<", SL_PP: "<<djack_t(SL_AP_QED[ind_QED]).ave_err()<<endl;
	  cout<<plots_path<<", D(Z)A/(Z)A: "<<djack_t(DZA_QED_rel[ind_QED]).ave_err()<<endl;
	  cout<<plots_path<<", D(Z)A: "<<djack_t(DZA_QED_rel[ind_QED]*ZA[ind_QCD]).ave_err()<<endl;
	}
      
      for(size_t iquark=0;iquark<4;iquark++)
	cout<<"Kritical kappa shift for iens "<<iens<<" "<<deltam_cr[iens].ave_err()<<" quark "<<iquark<<": "<<djack_t(deltam_cr[iens]*sqr(all_qpars[iquark].eq)*e2).ave_err()<<endl;
    }
}

/////////////////////////////////////////////////////////// dml_ren /////////////////////////////////////////////

namespace dml
{
  const size_t ncont_extrap=2;
  const size_t nchir_variations=2;
  const size_t nFSE_variations=2;
  
  index_t ind_syst({{"Input",ninput_an},
		    {"Chir",nchir_variations},
		    {"FSE",nFSE_variations},
		    {"Cont",ncont_extrap}});
  enum syst{c_input,c_chir,c_FSE,c_cont};
  template <syst comp> int case_of(int isyst){return ind_syst(isyst)[comp];}
}

template <class Tpars> Tpars FSE_dml_ren(const Tpars &C,const Tpars &L3dep,const double &L,const Tpars &powL)
{return C*L3dep*pow(L,powL);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dml_ren(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const Tpars &powL)
{
  Tpars xi=xi_fun(B0,aml,f0);
  return C*(1.0+Kpi*xi+Kpi*xi*xi+a*a*adep)+FSE_dml_ren(C,L3dep,L,powL);
}

//! perform the fit to the continuum limit of QED
dboot_t cont_chir_fit_dml_ren(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t isyst,bool cov_flag,const vector<string> &beta_list)
{
  using namespace dml;
  
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  const ave_err_t adep_ml_guess={0,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0,0.001};
  const ave_err_t powL_guess={-2,0.001};
  const string yaxis_title="$$\\delta m_l^{ren}";
  
  ave_err_t C_guess(0.0012,0.0012);
  ave_err_t KPi_guess(2.0,0.2);
  ave_err_t adep_guess={0,0.1};
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_QED",C_guess.ave(),C_guess.err());
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
  pars.iL4dep=boot_fit.add_fit_par(pars.L4dep,"powL",powL_guess.ave(),powL_guess.err());
  boot_fit.fix_par_to(pars.iL4dep,-3.0);
  
  boot_fit.fix_par_to(pars.iK2Pi,0.0);
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  if(case_of<c_FSE>(isyst)==0) boot_fit.fix_par_to(pars.iL3dep,0.0);
  //if(case_of<c_cont>(isyst)==0)
  //boot_fit.fix_par_to(pars.iadep,0.0);
  if(case_of<c_chir>(isyst)==1) boot_fit.fix_par_to(pars.iKPi,0.0);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[isyst](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_dml_ren(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],p[pars.iL4dep]);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_dml_ren(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,pars.L4dep);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  if(case_of<c_FSE>(isyst)==1)
    {
      grace::default_color_scheme={grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
      grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
    }
  else
    {
      grace::default_color_scheme={grace::RED,grace::BLUE,grace::GREEN4,grace::VIOLET};
      grace::default_symbol_scheme={grace::DIAMOND,grace::DIAMOND,grace::DIAMOND};
    }
  
  plot_chir_fit(path,ext_data,pars,
		[&pars]
		(double x,size_t ib)
		{return cont_chir_ansatz_dml_ren<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),
		     inf_vol,pars.L3dep.ave(),pars.L4dep.ave());},
		bind(cont_chir_ansatz_dml_ren<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,
		     inf_vol,pars.L3dep,pars.L4dep),
		[&ext_data,&pars,&B0,&f0]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_dml_ren(pars.C,pars.L3dep,ext_data[idata].L,pars.L4dep));},
		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  return phys_res;
}

//! compute the correction to the bare masses needed
index_t ind_adml;
dbvec_t adml_bare;
void compute_adml_bare()
{
  ind_adml.set_ranges({{"Input",ninput_an},{"Ens",nens_used}});
  adml_bare.resize(ind_adml.max());
  
  dbvec_t dml_ren_contlin(dml::ind_syst.max());
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      vector<cont_chir_fit_data_t> data_dml_all,data_dml_use_for_L;
      prepare_az(input_an_id);
      
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  size_t ind=ind_adml({input_an_id,iens});
	  ens_pars_t &ens=ens_pars[iens];
	  bi=jack_index[input_an_id][ens.iult];
	  
	  size_t ib=ens.ib;
	  dboot_t a=1/lat_par[input_an_id].ainv[ib];
	  dboot_t ZP=lat_par[input_an_id].Z[ib];
	  
	  size_t ind_ens_K=ind_ens_QCD_mes({iens,iK});
	  size_t ind_ens_Kplus=ind_ens_QED_mes({iens,iKPlus});
	  size_t ind_ens_K0=ind_ens_QED_mes({iens,iK0});
	  const double phys_dM2K=sqr(MKPLUS)-sqr(MK0);
	  
	  dboot_t QED_dM2K;
	  QED_dM2K=dboot_t(bi,SL_PP_QED[ind_ens_Kplus]-SL_PP_QED[ind_ens_K0]);
	  QED_dM2K*=e2*2*dboot_t(bi,M[ind_ens_K]);
	  QED_dM2K-=dboot_t(bi,FVE_M2(M[ind_ens_K],ens.L));
	  QED_dM2K/=sqr(a);
	  
	  dboot_t QCD_dM2K=phys_dM2K-QED_dM2K;
	  
	  //compute the proportionality factor between the mass slope of
	  //dM2K and aml_bare, needed to know by how much to multiply any mass correlator
	  dboot_t QCD_dM2K_over_adm;
	  QCD_dM2K_over_adm=dboot_t(bi,(SL_PP_MASS[ind_ens_Kplus]-SL_PP_MASS[ind_ens_K0])*2*M[ind_ens_K])/sqr(a);
	  adml_bare[ind]=QCD_dM2K/QCD_dM2K_over_adm;
	  
	  dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*ZP*(6.0*log(mu_MS*a)-22.596)/(32.0*sqr(M_PI)));
	  dboot_t ml=ens_pars[iens].aml/ZP/a;
	  dboot_t dml_ren=adml_bare[ind]/ZP/a-ml/Z_QED;
	  cout<<"dm_ren["<<ind<<"]: "<<dml_ren.ave_err()<<endl;
	  
	  dboot_t dum;
	  dum=0.0;
	  cont_chir_fit_data_t temp(ens.aml,ens.aml,dum,ib,ens.L,dml_ren,dml_ren);
	  data_dml_all.push_back(temp);
	  if(ens.use_for_L) data_dml_use_for_L.push_back(temp);
	}
      
      //extrap
      if(0)
      for(size_t ichir=0;ichir<dml::nchir_variations;ichir++)
	for(size_t iFSE=0;iFSE<dml::nFSE_variations;iFSE++)
	  for(size_t icont=0;icont<dml::ncont_extrap;icont++)
	    {
	      size_t isyst=dml::ind_syst({input_an_id,ichir,iFSE,icont});
	      cout<<dml::ind_syst.descr(isyst)<<endl;
	      
	      vector<cont_chir_fit_data_t> &data_dml=((iFSE==1)?data_dml_all:data_dml_use_for_L);
	      dml_ren_contlin[isyst]=cont_chir_fit_dml_ren(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,
							   data_dml,lat_par[input_an_id].ml,combine("plots_dml/cont_chir_an%zu.xmg",isyst),isyst,false,beta_list);
	    }
    }
}

///////////////////////////////////////////////////// nasty diagram ////////////////////////////////////////////

const size_t NPROJ=1;
const size_t nw=9;
const size_t norie=2; //number of orientation of the meson
const size_t nrev=2; //number of possible reversion
const size_t nqins=3; //number of quarks inserted: 0(none), 1 or 2
const size_t nprocess=6;
const size_t nrlep=2;
const size_t nproj=1;
index_t ind_hl_corr;

//! load hl correlations
djvec_t load_hl(size_t iproc,size_t iw,size_t iproj,const int *orie_par,/* const int *rev_par,*/size_t qins,const ens_pars_t &ens,const string &name)
{
  size_t T=ens.T;
  djvec_t out(T);
  out=0;
  
  size_t iQED_mes=iQED_mes_of_proc[iproc];
  size_t irev=QED_mes_pars[iQED_mes].irev;
  
  size_t n=0;
  //for(size_t irev=0;irev<nrev;irev++)
    for(size_t r2=0;r2<nr;r2++)
      for(size_t orie=0;orie<norie;orie++)
	for(size_t rl=0;rl<nr;rl++)
	  //for(size_t ri=0;ri<2;ri++)
  	    {
	      size_t ri=0;
	      size_t ic=ind_hl_corr({iproc,qins,irev,r2,orie,rl,iw,iproj,ri});
	      djvec_t corr=read_djvec(ens.path+"/data/corr_hl",T,ic);
	      
	      //insertion on
	      if(ri==RE) // not doing this to remove mizing! and r2==rl) //nb keeping r2 and rl identical
		{
		  //double r=rev_par[irev];
		  //if(qins==0) NOOOO dependency! we must keep the sign in place!!!!!
		  int r=1;
		  
		  out+=orie_par[orie]*r*corr;
		  n++;
		}
	      
	      if(0 and name!="")
		{
		  string path=ens.path+"/plots_hl/"+name+"_proc_"+to_string(iproc)+"_orie_"+to_string(orie)+"_qins_"+to_string(qins)+"_qrev_"+to_string(irev+1)+
		    "_r2_"+to_string(r2)+"_rl_"+to_string(rl)+"_ri_"+to_string(ri)+".xmg";
		  grace_file_t fout(path);
		  fout.write_vec_ave_err(corr.ave_err());
		  fout.set_title("iw="+to_string(iw)+" "+to_string(ic));
		}
	    }
    //cout<<n<<endl;
  return out.symmetrized(1)/n*pow(ens.L,3);
}

//! subtract the around-the-world effect
djvec_t hl_corr_subtract_around_world(const djvec_t &in,const djack_t &M)
{
  djvec_t out=in;
  djack_t fp=exp(M),fm=exp(-M);
  
  for(size_t t=1;t<out.size()-1;t++)
    out[t]=(in[t]+(in[t-1]*fp-in[t+1]*fm)/(fp-fm))/2;
  
  return out;
}

//! load the correlation and correct for around-the-world effect
valarray<djvec_t> load_and_correct_hl(size_t iproc,size_t iw,size_t iproj,const int *orie_par,/*const int *rev_par,*/size_t iens,const string &name)
{
  ens_pars_t &ens=ens_pars[iens];
  size_t iQCD_mes=QED_mes_pars[iQED_mes_of_proc[iproc]].iQCD;
  const djack_t M0=M[ind_ens_QCD_mes({iens,iQCD_mes})];
  djack_t mismatch=M0-ens.MMes[iQCD_mes];
  
  valarray<djvec_t> out(3);
  for(size_t qins=0;qins<3;qins++)
    {
      djvec_t precorr=load_hl(iproc,iw,iproj,orie_par,/*rev_par,*/qins,ens,name);
      djvec_t postsub=hl_corr_subtract_around_world(precorr,M0);
      
      djvec_t postmism=postsub;
      for(size_t t=0;t<postmism.size();t++) postmism[t]*=exp(mismatch*t);
      
      grace_file_t plot(ens.path+"/plots_hl/"+name+"_proc_"+to_string(iproc)+"_qins_"+to_string(qins)+".xmg");
      plot.set_subtitle("M0["+to_string(iens)+","+to_string(iproc)+"]: "+smart_print(M0.ave_err())+" input["+to_string(iQCD_mes)+"]: "+to_string(ens.MMes[iQCD_mes]));
      plot.write_vec_ave_err(precorr.ave_err());
      plot.set_legend("Raw");
      plot.write_vec_ave_err(postsub.ave_err());
      plot.set_legend("Corrected Around-World");
      plot.write_vec_ave_err(postmism.ave_err());
      plot.set_legend("Corrected for Mismatch of "+smart_print(mismatch.ave_err()));
      
      out[qins]=postmism;
    }
  
  return out;
}

//! load all correlation hl
vector<djvec_t> jLO_A_bare,jQED_V_bare,jQED_A_bare;
index_t ind_ens_proc;
void load_all_hl()
{
  ind_ens_proc.set_ranges({{"Ens",nens_used},{"Proc",nprocess}});
  size_t nens_proc=ind_ens_proc.max();
  jLO_A_bare.resize(nens_proc);
  jQED_V_bare.resize(nens_proc);
  jQED_A_bare.resize(nens_proc);
  
  const size_t nw=9;
  ind_hl_corr.set_ranges({{"NProcess",nprocess},{"NQIns",nqins},{"NRev",nrev},{"Nr2",nr},{"NOrie",norie},{"Nrlep",nrlep},{"NWeak",nw},{"NProj",nproj},{"RI",2}});
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      const size_t iqVi_lVi=0,iqV0_lV0=1; //all null
      const size_t iqAi_lVi=6,iqA0_lV0=7;
      //const int iqAi_lAi=2,iqA0_lA0=3,iqVi_lAi=4,iqV0_lA0=5; //redundant because we put 1-g5 in the lepton side
      const size_t ipV0=0;
      //const int unk[2]={+2,+0};
      //const int unk2[2]={0,+2};
      const int evn[2]={+1,+1};
      //const int odd[2]={+1,-1};
      for(size_t iproc=0;iproc<nprocess;iproc++)
	{
	  size_t iQED_mes=iQED_mes_of_proc[iproc];
	  valarray<djvec_t>
	    qVi_lVi_pV0_allins=load_and_correct_hl(iproc,iqVi_lVi,ipV0,evn,/*odd,*/iens,"ViVi"),
	    qV0_lV0_pV0_allins=load_and_correct_hl(iproc,iqV0_lV0,ipV0,evn,/*odd,*/iens,"V0V0"),
	    qAi_lVi_pV0_allins=load_and_correct_hl(iproc,iqAi_lVi,ipV0,evn,/*odd,*/iens,"AiVi"),
	    qA0_lV0_pV0_allins=load_and_correct_hl(iproc,iqA0_lV0,ipV0,evn,/*odd,*/iens,"A0V0"),
	    qV_lV_pV0_allins=qV0_lV0_pV0_allins+qVi_lVi_pV0_allins,
	    qA_lV_pV0_allins=qA0_lV0_pV0_allins+qAi_lVi_pV0_allins;
	  
	  size_t ind=ind_ens_proc({iens,iproc});
	  jLO_A_bare[ind]=-qA0_lV0_pV0_allins[0]; //minus because (V-A)*(V-A)=2*(VV-AV) but the two is inside
	  double eq1=QED_mes_pars[iQED_mes].eq1,eq2=QED_mes_pars[iQED_mes].eq2; //the same should be done here, but we do below
	  jQED_V_bare[ind]=eq1*qV_lV_pV0_allins[1]+eq2*qV_lV_pV0_allins[2];
	  jQED_A_bare[ind]=eq1*qA_lV_pV0_allins[1]+eq2*qA_lV_pV0_allins[2];
	}
    }
}

// void fse(double mlep,const dboot_t &MPS)
// {
//   double MW=80.385;
//   double pi2=M_PI*M_PI;
//   dboot_t rl=mlep/MPS,rl2=rl*rl;
//   dboot_t bIR=1/(8*pi2)*(1+rl2)/(1-rl2)*log(rl2);
//   dboot_t b0=1/(16*pi2)*(-log(MW/MPS)+37/12.0+log(rl2)*(2*(1-3*rl2)+(1+rl2)*log(rl2))/(2*(1+rl2))+(1+rl2)/(1-rl2)*log(rl2)*(gammaeul-log(4*M_PI))+2*(1+rl2)*(k31+k32));
  
// }

//! compute the correction to the process
void compute_corr(size_t iproc)
{
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      size_t ib=ens.ib;
      size_t ind=ind_ens_proc({iens,iproc});
      
      for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
	{
	  prepare_az(input_an_id);
	  bi=jack_index[input_an_id][ens.iult];
	  
	  dbvec_t LO=Zv[ib]*dbvec_t(bi,jLO_A_bare[ind]);
	  dbvec_t QED=-dbvec_t(bi,jQED_A_bare[ind])*Zv[ib]+dbvec_t(bi,jQED_V_bare[ind])*Za[ib]; //minus as explained before
	  LO.ave_err().write(combine("%s/plots_hl/LO_iproc%zu_ian%zu.xmg",ens.path.c_str(),iproc,input_an_id));
	  QED.ave_err().write(combine("%s/plots_hl/QED_iproc%zu_ian%zu.xmg",ens.path.c_str(),iproc,input_an_id));
	  
	  dbvec_t rat_ext=QED/LO;
	  rat_ext[rat_ext.size()-1]=rat_ext[0]=0.0; //set to zero the contact term
	  rat_ext.ave_err().write(combine("%s/plots_hl/QED_LO_ratio_iproc%zu_ian%zu.xmg",ens.path.c_str(),iproc,input_an_id));
	  
	  //DZA_QED_rel[ind_QED];
	  
	}
    }
}

int main(int narg,char **arg)
{
  int start=time(0);
  
  cout.precision(16);
  cout<<zeta(0.27138338825)<<endl;
  cout<<endl;
  cout<<zeta(0)<<endl;
  // exit(0);
  
  initialize(narg,arg);
  
  compute_basic_slopes();
  compute_adml_bare();
  load_all_hl();
  
  for(size_t iproc=0;iproc<nprocess;iproc++)
    compute_corr(iproc);
  
  cout<<endl<<"Total time: "<<time(0)-start<<" s"<<endl;
  
  return 0;
}
