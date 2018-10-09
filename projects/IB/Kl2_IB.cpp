#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

const bool EXCLUDE_SMALL_VOLS=false;
const bool EXCLUDE_HIGH_MASSES=false;

//#define XI

#include <set>

#include <tranalisi.hpp>
#include <Kl2_IB_FSE.hpp>
#include <common.hpp>

const size_t nfit_range_variations=2;
const int frange_var[2]={0,-1};

////////////////////////////////////////////////// quarks /////////////////////////////////////////////////

class qpars_t
{
public:
  size_t iq;
  double dm;
  double eq;
  qpars_t(size_t iq,double dm,double eq) : iq(iq),dm(dm),eq(eq) {}
};
const int iqU=0,iqD=0,iqS=1,iqC=2;
const int dmU=-1,dmD=+1,dmS=0,dmC=0;
const qpars_t qU(iqU,dmU,eu),qD(iqD,dmD,ed),qS(iqS,dmS,es),qC(iqC,dmC,ec);
const vector<qpars_t> all_qpars({qU,qD,qS,qC});

////////////////////////////////////////////////// QCD mesons //////////////////////////////////////////////

//! hold name and consittuent quark for a given meson
class QCD_mes_pars_t
{
public:
  size_t iq1; //!< index of quark 1
  size_t iq2; //!< index of quark 2
  size_t itint; //!< which interval range this meson has to use (0=Pi, 1=K, 2=D,Ds)
  string name; //!< name of the quark
  QCD_mes_pars_t(const qpars_t &q1,const qpars_t &q2,size_t itint,const string &name) : iq1(q1.iq),iq2(q2.iq),itint(itint),name(name) {}
};

//! for pure QCD meson we take neutral combinations
const vector<QCD_mes_pars_t> QCD_mes_pars({{qD,qD,0,"Pi"},{qD,qS,1,"K"},{qC,qD,2,"D"},{qC,qS,2,"Ds"}});
enum{iPi,iK,iD};
const size_t &nQCD_mes=QCD_mes_pars.size();

/////////////////////////////////////////////////// QED mesons /////////////////////////////////////

//! hold name and consittuent quark for a given QED meson
class QED_mes_pars_t
{
public:
  const size_t iq1; //!< index of quark 1
  const size_t iq2; //!< index of quark 2
  const double dm1; //!< contribution by which scalar insertion on quark 1 must be added
  const double dm2; //!< contribution by which scalar insertion on quark 2 must be added
  const double eq1; //!< charge of quark 1
  const double eq2; //!< charge of quark 2
  const size_t irev; //!< quark line to be reversed
  const size_t iQCD; //!< index of pure QCD meson corresponding to this meson
  const double M; //!< mass of the meson
  const string name; //!< name of the meson
  QED_mes_pars_t(const qpars_t &q1,const qpars_t &q2,const size_t irev,const size_t iQCD,const double M,const string &name) :
    iq1(q1.iq),iq2(q2.iq),dm1(q1.dm),dm2(q2.dm),eq1(q1.eq),eq2(q2.eq),irev(irev),iQCD(iQCD),M(M),name(name) {}
};

const size_t QREV1=0; //!< value to revert quark 1
const size_t QREV2=1; //!< value to revert quark 2
const vector<QED_mes_pars_t> QED_mes_pars({{qD,qU,QREV2,0,MPPLUS,"PiMinus"},
                                           {qU,qS,QREV1,1,MKPLUS,"KMinus"},
					   {qD,qS,QREV1,1,MK0,"K0bar"},
					   {qD,qC,QREV2,2,MDPLUS,"DMinus"},
					   {qU,qC,QREV2,2,MD0,"D0bar"},
					   {qS,qC,QREV2,3,MDS,"DsMinus"}});
enum{iPiMinus,iKMinus,iK0bar,iDMinus,iD0bar,iDsMinus};
const size_t &nQED_mes=QED_mes_pars.size();

////////////////////////////////////////////////////////////////////

const size_t nleps=2; //!< number of leptons
const size_t nmes_tint=3; //!< number of time intervals
const size_t nqmass=3; //!< number of quark mass

class ens_pars_t
{
public:
  size_t iult; //!< input in the ultimate file
  size_t ib; //!< beta index
  size_t T; //!< time extent
  size_t L; //!< spatial size
  double am[nqmass]; //! quark mass
  double &aml=am[0]; //!< bare light quark mass
  double &ams=am[1]; //!< bare strange quark mass
  double &amc=am[2]; //!< bare charm quark mass
  double aMLep[nleps]; //!< mass of leptons
  double aMMes[4]; //!< mass of mesons (0=Pi, 1=K, 2=D, 3=Ds)
  bool use_for_L; //!< use for FSE analysis
  bool use_for_chir; //!< use for chiral analysis
  //bool use_for_a; //!< use for cont analysis
  string path; //!< path (name)
  
  vector<size_t> tmin,tmax; //!< range of fit
  ens_pars_t() : tmin(nmes_tint),tmax(nmes_tint) {}
};
vector<ens_pars_t> ens_pars; //!< parameters of all ensemble
size_t nens_used; //!< number of ensemble used

const size_t nr=2; //!< number of r
const index_t ind_2pts({{"NMass",nqmass},{"NMass",nqmass},{"Nr",nr},{"RI",2}});

/* Processes are
   Pi->Mu
   K->Mu
   D->Mu
   D->Tau
   Ds->Mu
   Ds->Tau */

//! The decaying meson is always the negatively charged one
vector<size_t> iQED_mes_of_proc({0,1,3,3,5,5}); //!< index of the QED meson corresponding to a given process from the QED_mes_pars
const double MLep[2]={0.1056583745,1.77682};
vector<size_t> iMLep_of_proc({0,0,0,1,0,1}); //!< index of the lepton mass corresponding to a given process

//! read or write z0 for FSE
vector<double> z0;
void prepare_z0()
{
  string path_z0("z0.dat");
  if(not file_exists(path_z0))
    {
      cout<<"Preparing z0"<<endl;
      raw_file_t(path_z0,"w").bin_write(z0=zeta_FSE(0.0));
    }
  else
    {
      cout<<"Reading z0"<<endl;
      z0.resize(nZ_FSE);
      raw_file_t(path_z0,"r").bin_read(z0);
    }
}

//! initialize everything
void initialize(int narg,char **arg)
{
  prepare_z0();
  
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
  
  input.expect({"Ens","beta","aml","ams","amc","L","UseForL","T","TPi","TK","TD","path","aMLep0","aMLep1","aMMes0","aMMes1","aMMes2","aMMes3"});
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      
      input.read(ens.iult);
      input.read(ens.ib);
      input.read(ens.aml);
      input.read(ens.ams);
      input.read(ens.amc);
      input.read(ens.L);
      input.read(ens.use_for_L);
      //input.read(ens.use_for_a);
      //ens.use_for_a=(ens.ib!=0);
      input.read(ens.T);
      for(size_t itint=0;itint<nmes_tint;itint++)
	{
	  input.read(ens.tmin[itint]);
	  input.read(ens.tmax[itint]);
	}
      input.read(ens.path);
      for(size_t ilep=0;ilep<nleps;ilep++) input.read(ens.aMLep[ilep]);
      for(size_t iQCD_mes=0;iQCD_mes<nQCD_mes;iQCD_mes++) input.read(ens.aMMes[iQCD_mes]);
      
      //replace use_for_L with a made-up one
      if(EXCLUDE_SMALL_VOLS) ens.use_for_L=(ens.L*ens.aMMes[0]>4.5);
      else                   ens.use_for_L=true;
      //switch chir inclusion
      if(EXCLUDE_HIGH_MASSES) ens.use_for_chir=(ens.aMMes[0]*lat_par[0].ainv[ens.ib].ave()<=0.350);
      else                    ens.use_for_chir=true;
    }
}

//! read a single vector, for a specific mass and r, real or imaginary
djvec_t read(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,size_t r,size_t reim)
{
  string path=combine("%s/data/corr%s",ens.path.c_str(),what);
  djvec_t out=read_djvec(path,ens.T,ind_2pts({iq1,iq2,r,reim}));
  return out;
}

//! read a combination of r and return appropriately symmetrized in time and r
djvec_t read(const char *what,const ens_pars_t &ens,int tpar,size_t iq1,size_t iq2,int rpar,size_t reim)
{
  djvec_t o(ens.T);
  o=0.0;
  
  for(size_t r=0;r<nr;r++) o+=read(what,ens,iq1,iq2,r,reim)*((r==0)?1:rpar);
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
  const string ens_qpath=ens.path+"/plots_mcrit";
  
  const djvec_t V0P5_LL=read_VP("LL",ens,iq,iq,-1,IM);
  const djvec_t V0P5_0M=read_VP("0M",ens,iq,iq,-1,IM);
  const djvec_t V0P5_0T=read_VP("0T",ens,iq,iq,-1,IM);
  const djvec_t num_deltam_cr=forward_derivative(djvec_t(V0P5_LL+2.0*djvec_t(V0P5_0M+V0P5_0T)));
  num_deltam_cr.ave_err().write(combine("%s/num_deltam_cr.xmg",ens_qpath.c_str()));
  
  const djvec_t V0P5_0P=read_VP("0P",ens,iq,iq,+1,RE);
  const djvec_t den_deltam_cr=forward_derivative(V0P5_0P);
  den_deltam_cr.ave_err().write(combine("%s/den_deltam_cr.xmg",ens_qpath.c_str()));
  
  const size_t itint=QCD_mes_pars[iQCD_mes].itint;
  const djack_t deltam_cr=constant_fit(djvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),ens.tmin[itint],ens.tmax[itint],combine("%s/deltam_cr_t.xmg",ens_qpath.c_str()));
  
  return deltam_cr;
}

//! read QED corrections (inner diagrams)
djvec_t read_QED(const ens_pars_t &ens,size_t iQED_mes,const djack_t &deltam_cr,const djvec_t &c_LO,
		 djvec_t(*read)(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim),const char *name)
{
  const QED_mes_pars_t &pars=QED_mes_pars[iQED_mes];
  string ens_qpath=ens.path+"/plots_"+pars.name;
  
  size_t iq1=pars.iq1;
  size_t iq2=pars.iq2;
  double eq1=pars.eq1;
  double eq2=pars.eq2;
  double am1=ens.am[iq1];
  double am2=ens.am[iq2];
  const double a=1/lat_par[0].ainv[ens.ib].ave();
  
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
  djvec_t c=c_LL+c_0T1+c_0M1+c_0T2+c_0M2;
  
  djvec_t c_0P1=-(read("0P",ens,iq2,iq1,-1,IM)*eq1*eq1); //remind that i missing means
  djvec_t c_0P2=-(read("0P",ens,iq1,iq2,-1,IM)*eq2*eq2); //to take imag part changed of sign
  djvec_t(c_0P1/c_LO).ave_err().write(combine("%s/%s_0P1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0P2/c_LO).ave_err().write(combine("%s/%s_0P2.xmg",ens_qpath.c_str(),name));
  djvec_t d=deltam_cr*(c_0P1+c_0P2);
  
  //subtract the bare quark mass eq.85 of PRD 2013
  double Z_QED=(6.0*log(mu_MS*a)-22.5954)/(16.0*sqr(M_PI));

  djack_t ZP_fact;
  ZP_fact.fill_gauss({1.5,0.2,23492});
  
  //the insertion is of the scalar density, but
  //the correction is given by the lagrangian
  //insertion, where -S is present
  djvec_t c_0S1=-read("0S",ens,iq2,iq1,1,RE)*am1*sqr(eq1)*Z_QED*ZP_fact;
  djvec_t c_0S2=-read("0S",ens,iq1,iq2,1,RE)*am2*sqr(eq2)*Z_QED*ZP_fact;
  
  djvec_t(c_0S1/c_LO).ave_err().write(combine("%s/%s_0S1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0S2/c_LO).ave_err().write(combine("%s/%s_0S2.xmg",ens_qpath.c_str(),name));
  djvec_t e=c_0S1+c_0S2;
  
  return c+d+e;
}

//! read MASS corrections
djvec_t read_MASS(const ens_pars_t &ens,size_t iQED_mes,const djvec_t &c_LO,
		  djvec_t(*read)(const char *what,const ens_pars_t &ens,size_t iq1,size_t iq2,int rpar,size_t reim),const char *name)
{
  string ens_qpath=ens.path+"/plots_"+QED_mes_pars[iQED_mes].name;
  const QED_mes_pars_t &pars=QED_mes_pars[iQED_mes];
  size_t iq1=pars.iq1;
  size_t iq2=pars.iq2;
  //the insertion is of the scalar density, but
  //the correction is given by the lagrangian
  //insertion, where -S is present
  djvec_t c_0S1=-read("0S",ens,iq2,iq1,1,RE)*pars.dm1;
  djvec_t c_0S2=-read("0S",ens,iq1,iq2,1,RE)*pars.dm2;
  
  djvec_t(c_0S1/c_LO).ave_err().write(combine("%s/%s_0S1.xmg",ens_qpath.c_str(),name));
  djvec_t(c_0S2/c_LO).ave_err().write(combine("%s/%s_0S2.xmg",ens_qpath.c_str(),name));
  
  return c_0S1+c_0S2;
}

//slopes and Z for all mesons and masses
djvec_t jZP,jZA;
djvec_t jaM,jxi;
djvec_t jDZA_QED_rel,jDZA_MASS_rel;
djvec_t jDM_QED,jDM_MASS;

index_t ind_ens_QCD_mes;
index_t ind_ens_QED_mes;
djvec_t deltam_cr;
vector<djvec_t> jAP_LO_exp_removed; //!< used to test circle
void compute_basic_slopes()
{
  ind_ens_QCD_mes.set_ranges({{"Ens",nens_used},{"Mes",nQCD_mes}});
  ind_ens_QED_mes.set_ranges({{"Ens",nens_used},{"Mes",nQED_mes}});
  size_t nens_QCD_mes=ind_ens_QCD_mes.max();
  size_t nens_QED_mes=ind_ens_QED_mes.max();
  
  deltam_cr.resize(nens_used);
  jZP.resize(nens_QCD_mes);
  jZA.resize(nens_QCD_mes);
  jaM.resize(nens_QCD_mes);
  jxi.resize(nens_QCD_mes);
  //
  jDZA_QED_rel.resize(nens_QED_mes);
  jDZA_MASS_rel.resize(nens_QED_mes);
  //
  jDM_QED.resize(nens_QED_mes);
  jDM_MASS.resize(nens_QED_mes);
  //
  
  vector<djvec_t> jPP_LO(nens_QCD_mes);
  vector<djvec_t> jPP_MASS(nens_QED_mes);
  vector<djvec_t> jPP_QED(nens_QED_mes);
  //
  vector<djvec_t> jAP_LO(nens_QCD_mes);
  vector<djvec_t> jAP_QED(nens_QED_mes);
  vector<djvec_t> jAP_MASS(nens_QED_mes);
  jAP_LO_exp_removed.resize(nens_QCD_mes);

  FILE* pFile;

  pFile=fopen("jDM_MASS_K0.dat","w+");
  fclose(pFile);
  
  //load everything
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      deltam_cr[iens]=compute_deltam_cr(ens,ilight,iPi);
      size_t T=ens.T,TH=T/2;
      
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
	  //load LO for PP
	  jPP_LO[ind_QCD]=read_PP("00",ens,iq1,iq2,1,RE);
	  cout<<"jPP: "<<jPP_LO[ind_QCD][1]<<endl; //check normalization
	  //load corrections for PP
	  jPP_MASS[ind_QED]=read_MASS(ens,iQED_mes,jPP_LO[ind_QCD],read_PP,"PP");
	  jPP_QED[ind_QED]=read_QED(ens,iQED_mes,deltam_cr[iens],jPP_LO[ind_QCD],read_PP,"PP");
	  //load LO for AP
	  jAP_LO[ind_QCD]=read_AP("00",ens,iq1,iq2,1,RE);
	  jAP_LO[ind_QCD].ave_err().write(plots_path+"/corr_AP_LO.xmg");
	  //load corrections for AP
	  jAP_MASS[ind_QED]=read_MASS(ens,iQED_mes,jAP_LO[ind_QCD],read_AP,"AP");
	  jAP_QED[ind_QED]=read_QED(ens,iQED_mes,deltam_cr[iens],jAP_LO[ind_QCD],read_AP,"AP");
	  
	  //fit separately M in each channel, to be improved
	  djack_t ZAP,ZPP;
	  djack_t SL_AP_QED,SL_PP_QED;
	  djack_t DZ_AP_QED_rel,DZ_PP_QED_rel;
	  djack_t SL_AP_MASS,SL_PP_MASS;
	  djack_t DZ_AP_MASS_rel,DZ_PP_MASS_rel;
	  two_pts_with_ins_ratio_fit(ZAP         ,jaM[ind_QCD],DZ_AP_QED_rel,SL_AP_QED,jAP_LO[ind_QCD],jAP_QED[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_AP_LO1.xmg",plots_path+"/slope_AP_QED.xmg",-1);
	  two_pts_with_ins_ratio_fit(ZAP         ,jaM[ind_QCD],DZ_AP_MASS_rel,SL_AP_MASS,jAP_LO[ind_QCD],jAP_MASS[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_AP_LO2.xmg",plots_path+"/slope_AP_MASS.xmg",-1);
	  //
	  two_pts_with_ins_ratio_fit(ZPP,jaM[ind_QCD],DZ_PP_QED_rel,SL_PP_QED,jPP_LO[ind_QCD],jPP_QED[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_PP_LO1.xmg",plots_path+"/slope_PP_QED.xmg");
	  two_pts_with_ins_ratio_fit(ZPP,jaM[ind_QCD],DZ_PP_MASS_rel,SL_PP_MASS,jPP_LO[ind_QCD],jPP_MASS[ind_QED],TH,tmin,tmax,
				     plots_path+"/effmass_PP_LO2.xmg",plots_path+"/slope_PP_MASS.xmg");
	  //
	  jZP[ind_QCD]=sqrt(ZPP);
	  jZA[ind_QCD]=ZAP/jZP[ind_QCD];
	  djack_t jaf=jZP[ind_QCD]*(ens.am[iq1]+ens.am[iq2])/sqr(jaM[ind_QCD]);
	  jxi[ind_QCD]=sqr(djack_t(jaM[ind_QCD]/(4*M_PI*jaf)));
	  cout<<"xi"<<iQCD_mes<<" "<<ens.path<<" "<<jxi[ind_QCD]<<endl;
	  //
	  djack_t DZ_P_MASS_rel=DZ_PP_MASS_rel/2.0;
	  djack_t DZ_P_QED_rel=DZ_PP_QED_rel/2.0;
	  //
	  jDZA_MASS_rel[ind_QED]=DZ_AP_MASS_rel-DZ_P_MASS_rel;
	  jDZA_QED_rel[ind_QED]=DZ_AP_QED_rel-DZ_P_QED_rel;
	  jDM_QED[ind_QED]=-SL_PP_QED; //change signed due to slope definition
	  jDM_MASS[ind_QED]=-SL_PP_MASS; //change signed due to slope definition
	  
	  cout<<plots_path<<", corr[0]: "<<jPP_LO[ind_QCD][0]<<endl;
	  cout<<plots_path<<", ZA: "<<jZA[ind_QCD].ave_err()<<endl;
	  cout<<plots_path<<", M_AP: "<<jaM[ind_QCD].ave_err()<<endl;
	  cout<<plots_path<<", SL_AP_QED: "<<djack_t(SL_AP_QED).ave_err()<<endl;
	  cout<<plots_path<<", SL_PP_MASS: "<<djack_t(SL_PP_MASS).ave_err()<<endl;
	  cout<<plots_path<<", (DZA/ZA)_QED: "<<djack_t(jDZA_QED_rel[ind_QED]).ave_err()<<endl;
	  cout<<plots_path<<", (DZA/ZA)_MASS: "<<djack_t(jDZA_MASS_rel[ind_QED]).ave_err()<<endl;
	  cout<<plots_path<<", DZ_AP_MASS_rel: "<<DZ_AP_MASS_rel.ave_err()<<endl;
	  cout<<plots_path<<", DZ_P_MASS_rel: "<<DZ_P_MASS_rel.ave_err()<<endl;
	  
	  //store for later test
	  jAP_LO_exp_removed[ind_QCD]=jAP_LO[ind_QCD];
	  for(int t=0;t<=(int)TH;t++)
	    jAP_LO_exp_removed[ind_QCD][t]/=exp(-t*jaM[ind_QCD])+exp(-((int)T-t)*jaM[ind_QCD]);
	}

      ////test////
      size_t ind_ens_K0bar=ind_ens_QED_mes({iens,iK0bar});

      for(size_t i=0;i<jDM_MASS[ind_ens_K0bar].size();i++)
	{
	  pFile=fopen("jDM_MASS_K0.dat","a+");
	  fprintf(pFile,"%lg\n",jDM_MASS[ind_ens_K0bar][i]);
	  fclose(pFile);
	}
      
      //print some info useful for retuning kappa for each quark
      for(size_t iquark=0;iquark<4;iquark++)
	cout<<"Kritical kappa shift for iens "<<iens<<" "<<deltam_cr[iens].ave_err()<<" quark "<<iquark<<": "<<djack_t(deltam_cr[iens]*sqr(all_qpars[iquark].eq)*e2).ave_err()<<endl;
    }
}

/////////////////////////////////////////////////////////// dml_ren /////////////////////////////////////////////

/* The logic used to fix dml is different from that used in the quark
mass paper.  Here we impose that each ensemble reproduce the physical
difference between the mass of charged and neutral K, by first
subtracting QED effects, than determining the coefficient needed to
make the QCD correction correspond to the residual. Finally, the QED
contribution arising from the original bare quark mass is subtracted
perturbatively.*/

//! holds the systematics for the estimate of dml_re
namespace dml
{
  const size_t ncont_extrap=2;
  const size_t nchir_variations=2;
  const size_t nFSE_variations=2;
  
  index_t ind_syst(
       {{"Input",ninput_an},
	{"Chir",nchir_variations},
	{"FSE",nFSE_variations},
	{"Cont",ncont_extrap}});
  enum syst{c_input,c_chir,c_FSE,c_cont};
  
  template <syst comp>
  int case_of(int isyst)
  {return ind_syst(isyst)[comp];}
}

//! finit size effects on delta ml
template <class Tpars>
Tpars FSE_dml_ren(const Tpars &C,const Tpars &L3dep,const double &L,const Tpars &powL)
{return C*L3dep*pow(L,powL);}

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dml_ren(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &aml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const Tpars &powL)
{
  Tpars xi=xi_fun(B0,aml,aml,f0);
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
  
  const ave_err_t adep_ml_guess={0.0,0.001};
  const ave_err_t K2Pi_guess={-0.001,0.001};
  const ave_err_t L3dep_guess={0.0,0.001};
  const ave_err_t powL_guess={-2.0,0.001};
  const string yaxis_title="$$\\delta m_l^{ren}";
  
  ave_err_t C_guess(0.0012,0.0012);
  ave_err_t KPi_guess(2.0,0.2);
  ave_err_t adep_guess={0.0,0.1};
  
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
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double dum,double ac,double L)
			 {return cont_chir_ansatz_dml_ren(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],p[pars.iL4dep]);},cov_flag);
  
  double a_cont=0.0;
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
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_dml_ren(pars.C,pars.L3dep,ext_data[idata].L,pars.L4dep));},
		ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  return phys_res;
}

//! compute the correction to the bare masses needed
index_t ind_adml;
dbvec_t adml_bare;
dbvec_t MK;
void compute_adml_bare()
{
  ind_adml.set_ranges({{"Input",ninput_an},{"Ens",nens_used}});
  adml_bare.resize(ind_adml.max());
  MK.resize(ind_adml.max());

  FILE* pFile;

  pFile=fopen("adml_bare.dat","w+");
  fclose(pFile);

  FILE* ppFile;

  ppFile=fopen("MK.dat","w+");
  fclose(ppFile);
  
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
	  size_t ind_ens_Kminus=ind_ens_QED_mes({iens,iKMinus});
	  size_t ind_ens_K0bar=ind_ens_QED_mes({iens,iK0bar});
	  const double phys_dM2K=sqr(MKPLUS)-sqr(MK0);
	  
	  //compute the QED contribution
	  dboot_t QED_dM2K;
	  QED_dM2K=dboot_t(bi,jDM_QED[ind_ens_Kminus]-jDM_QED[ind_ens_K0bar]);
	  QED_dM2K*=e2*2*dboot_t(bi,jaM[ind_ens_K]);
	  QED_dM2K-=dboot_t(bi,FVE_M2(jaM[ind_ens_K],ens.L));
	  QED_dM2K/=sqr(a);
	  
	  //! subtract the QED contrubution
	  dboot_t QCD_dM2K=phys_dM2K-QED_dM2K;
	  
	  //compute the proportionality factor between the mass slope of
	  //dM2K and aml_bare, needed to know by how much to multiply any mass correlator
	  dboot_t QCD_dM2K_over_adm;
	  QCD_dM2K_over_adm=dboot_t(bi,(jDM_MASS[ind_ens_Kminus]-jDM_MASS[ind_ens_K0bar])*2*jaM[ind_ens_K])/sqr(a);
	  adml_bare[ind]=QCD_dM2K/QCD_dM2K_over_adm;

	  for(size_t iboot=0;iboot<adml_bare[ind].size();iboot++)
	    {
	      pFile=fopen("adml_bare.dat","a+");
	      fprintf(pFile,"%lg\n",adml_bare[ind][iboot]);
	      fclose(pFile);
	    }

	  MK[ind]=dboot_t(bi,jaM[ind_ens_K]);

	  for(size_t iboot=0;iboot<MK[ind].size();iboot++)
	    {
	      ppFile=fopen("MK.dat","a+");
	      fprintf(ppFile,"%lg\n",MK[ind][iboot]);
	      fclose(ppFile);
	    }
	  
	  //subtract the bare quark mass eq.85 of PRD 2013
	  dboot_t dml_ren=adml_bare[ind]/ZP/a;
	  cout<<"dm_ren["<<ind<<"]: "<<dml_ren.ave_err()<<endl;
	  
	  //! add to the fit
	  dboot_t dum;
	  dum=0.0;
	  cont_chir_fit_data_t temp(ens.aml,ens.ams,dum,ib,ens.L,dml_ren,dml_ren);
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
	      cout<<"dml (cont): "<<dml_ren_contlin[isyst].ave_err()<<endl;
	    }
    }
}

///////////////////////////////////////////////////// nasty diagram ////////////////////////////////////////////

//const size_t nw=9; //!< number of weak insertions
const size_t norie=2; //!< number of orientation of the meson
const size_t nrev=2; //!< number of possible reversion
const size_t nqins=3; //!< number of quarks inserted: 0(none), 1 or 2
const size_t nprocess_max=6; //!< number of processes computed
const size_t nprocess=2; //!< number of process to analyse
const size_t nrlep=2; //!< number of r for leptons
const size_t nproj=1; //!<number of projectors: 1, V0 only
index_t ind_hl_corr;

enum{STUDY_PI,STUDY_K_M_PI,STUDY_K};

//! holds the systematics for the estimate of hadroleptonic
namespace hl
{
  //! chiral extrapolation variations
  namespace chir
  {
    enum t{QUADRATIC,QUADRATICLOG};
    const vector<int> variations={QUADRATIC,QUADRATICLOG};
    const vector<string> tag={"QUADRATIC","QUADRATICLOG"};
    const size_t nvariations=variations.size();
  }

  namespace FSE
  {
    //! variations of data subtraction of FSE
    namespace SUB_STDEP
    {
      enum t{NO,YES};
      const vector<string> tag={"SUB_STDEP_NO","SUB_STDEP_YES"};
    }
    
    //! variations of fit of FSE
    namespace FIT_STDEP
    {
      enum t{NO,YES};
      const vector<string> tag={"FIT_STDEP_NO","FIT_STDEP_YES"};
    }
    
    //! consider two variations for the fit
    const size_t nvariations=2;
  }
  
  //! continuum extrapolation variations
  namespace cont
  {
    enum t{LINEAR,CONSTANT};
    const vector<int> variations={LINEAR,CONSTANT};
    const vector<string> tag={"LINEAR","CONSTANT"};
    const size_t nvariations=variations.size();
  }
  
  //! index of systematics
  index_t ind_syst({
      {"Input",ninput_an},
      {"Fit Range",nfit_range_variations},
      {"Chir",chir::nvariations},
      {"FSE",FSE::nvariations},
      {"Cont",cont::nvariations}});
  enum syst{c_input,c_frange,c_chir,c_FSE,c_cont};
  
  //! return the case of each systematic
  template <syst comp>
  int case_of(int isyst)
  {return ind_syst(isyst)[comp];}
  
  const vector<size_t> FSE_max_orders={1,2};
  index_t ind_an_ens_FSEmax_frange;
}

//! read hl correlations
djvec_t read_hl(size_t iproc,size_t iw,size_t iproj,const int *orie_par,size_t qins,const ens_pars_t &ens,const string &name,const array<int,2> &r2_weight={1,1},const array<int,2> rl_weight={1,1})
{
  size_t T=ens.T;
  djvec_t out(T);
  out=0;
  
  size_t iQED_mes=iQED_mes_of_proc[iproc];
  size_t irev=QED_mes_pars[iQED_mes].irev;
  
  //average the two r, the two orientations and the two r of the leptons
  size_t n=0;
  for(size_t r2=0;r2<nr;r2++)
    for(size_t orie=0;orie<norie;orie++)
      for(size_t rl=0;rl<nr;rl++)
	if(r2_weight[r2]*rl_weight[rl])
	  {
	  size_t ri=0; //only real part
	  size_t ic=ind_hl_corr({iproc,qins,irev,r2,orie,rl,iw,iproj,ri});
	  djvec_t corr=read_djvec(ens.path+"/data/corr_hl",T,ic);
	  
	  out+=orie_par[orie]*corr*r2_weight[r2]*rl_weight[rl];
	  n+=abs(r2_weight[r2]*rl_weight[rl]);
	  
	  // if(name!="")
	  //   {
	  //     string path=ens.path+"/plots_hl/"+name+"_proc_"+to_string(iproc)+"_orie_"+to_string(orie)+"_qins_"+to_string(qins)+"_qrev_"+to_string(irev+1)+
	  // 	"_r2_"+to_string(r2)+"_rl_"+to_string(rl)+"_ri_"+to_string(ri)+".xmg";
	  //     grace_file_t fout(path);
	  //     fout.write_vec_ave_err(corr.ave_err());
	  //     fout.set_title("iw="+to_string(iw)+" "+to_string(ic));
	  //   }
	}
  
  //if insertion is made, change sign, due to wrong sign in the program
  if(qins) out*=-1.0;
  
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
valarray<djvec_t> load_and_correct_hl(size_t iproc,size_t iw,size_t iproj,const int *orie_par,size_t iens,const string &name,const array<int,2> &r2_weight={1,1},const array<int,2> rl_weight={1,1})
{
  ens_pars_t &ens=ens_pars[iens];
  size_t iQCD_mes=QED_mes_pars[iQED_mes_of_proc[iproc]].iQCD;
  const djack_t M0=jaM[ind_ens_QCD_mes({iens,iQCD_mes})];
  djack_t mismatch=M0-ens.aMMes[iQCD_mes];
  
  //! load the three quark insertions
  valarray<djvec_t> out(3);
  for(size_t qins=0;qins<3;qins++)
    {
      //load and remove around the world
      djvec_t precorr=read_hl(iproc,iw,iproj,orie_par,qins,ens,name,r2_weight,rl_weight);
      cout<<"precorr["<<qins<<"]: "<<precorr[1]<<endl; //check normalization
      djvec_t postsub=hl_corr_subtract_around_world(precorr,M0);
      
      //remove mismatch in mass
      djvec_t postmism=postsub;
      for(size_t t=0;t<postmism.size();t++) postmism[t]*=exp(mismatch*t);
      
      //plot the comparison
      grace_file_t plot(ens.path+"/plots_hl/"+name+"_proc_"+to_string(iproc)+"_qins_"+to_string(qins)+".xmg");
      plot.set_subtitle("M0["+to_string(iens)+","+to_string(iproc)+"]: "+smart_print(M0.ave_err())+" input["+to_string(iQCD_mes)+"]: "+to_string(ens.aMMes[iQCD_mes]));
      plot.write_vec_ave_err(precorr.ave_err());
      plot.set_legend("Raw");
      plot.write_vec_ave_err(postsub.ave_err());
      plot.set_legend("Corrected Around-World");
      plot.write_vec_ave_err(postmism.ave_err());
      plot.set_legend("Corrected for Mismatch of "+smart_print(mismatch.ave_err()));
      
      //store the correlator after all corrections
      out[qins]=postmism;
    }
  
  return out;
}

//! load all correlation hl
vector<djvec_t> jLO_A_bare,jQED_V_bare,jQED_A_bare,jQED_qA0_lV0_bare,jQED_qV0_lV0_bare,jQED_qAi_lVi_bare,jQED_qVi_lVi_bare;
index_t ind_ens_proc;
void load_all_hl(const array<int,2> &r2_weight={1,1},const array<int,2> rl_weight={1,1})
{
  ind_ens_proc.set_ranges({{"Ens",nens_used},{"Proc",nprocess}});
  hl::ind_an_ens_FSEmax_frange.set_ranges(
			      {{"Input",ninput_an},
			       {"Ens",nens_used},
			       {"FSE max order",hl::FSE_max_orders.size()},
			       {"Fit range",nfit_range_variations}});
  size_t nens_proc=ind_ens_proc.max();
  jLO_A_bare.resize(nens_proc);
  jQED_V_bare.resize(nens_proc);
  jQED_A_bare.resize(nens_proc);
  /////////////////test Z_fact//////////////
  jQED_qA0_lV0_bare.resize(nens_proc);
  jQED_qV0_lV0_bare.resize(nens_proc);
  jQED_qAi_lVi_bare.resize(nens_proc);
  jQED_qVi_lVi_bare.resize(nens_proc);
  
  const size_t nw=9;
  ind_hl_corr.set_ranges({{"NProcess",nprocess_max},{"NQIns",nqins},{"NRev",nrev},{"Nr2",nr},{"NOrie",norie},{"Nrlep",nrlep},{"NWeak",nw},{"NProj",nproj},{"RI",2}});
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      const size_t iqVi_lVi=0,iqV0_lV0=1;
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
	    qVi_lVi_pV0_allins=load_and_correct_hl(iproc,iqVi_lVi,ipV0,evn,iens,"ViVi",r2_weight,rl_weight),
	    qV0_lV0_pV0_allins=load_and_correct_hl(iproc,iqV0_lV0,ipV0,evn,iens,"V0V0",r2_weight,rl_weight),
	    qAi_lVi_pV0_allins=load_and_correct_hl(iproc,iqAi_lVi,ipV0,evn,iens,"AiVi",r2_weight,rl_weight),
	    qA0_lV0_pV0_allins=load_and_correct_hl(iproc,iqA0_lV0,ipV0,evn,iens,"A0V0",r2_weight,rl_weight),
	    qV_lV_pV0_allins=qV0_lV0_pV0_allins+qVi_lVi_pV0_allins,
	    qA_lV_pV0_allins=qA0_lV0_pV0_allins+qAi_lVi_pV0_allins;
	  
	  size_t ind=ind_ens_proc({iens,iproc});
	  //leading order
	  jLO_A_bare[ind]=qA0_lV0_pV0_allins[0];
	  //include the photon attached to the quarks
	  double el=-1; //charge of electron and muon
	  double eq1=QED_mes_pars[iQED_mes].eq1,eq2=QED_mes_pars[iQED_mes].eq2;
	  jQED_V_bare[ind]=el*(eq1*qV_lV_pV0_allins[1]+eq2*qV_lV_pV0_allins[2]);
	  jQED_A_bare[ind]=el*(eq1*qA_lV_pV0_allins[1]+eq2*qA_lV_pV0_allins[2]);
	  /////////////////test Z_fact//////////////
	  jQED_qA0_lV0_bare[ind]=el*(eq1*qA0_lV0_pV0_allins[1]+eq2*qA0_lV0_pV0_allins[2]);
	  jQED_qV0_lV0_bare[ind]=el*(eq1*qV0_lV0_pV0_allins[1]+eq2*qV0_lV0_pV0_allins[2]);
	  jQED_qAi_lVi_bare[ind]=el*(eq1*qAi_lVi_pV0_allins[1]+eq2*qAi_lVi_pV0_allins[2]);
	  jQED_qVi_lVi_bare[ind]=el*(eq1*qVi_lVi_pV0_allins[1]+eq2*qVi_lVi_pV0_allins[2]);
	}
    }
}

// void fse(double MLep,const dboot_t &MPS)
// {
//   double MW=80.385;
//   double pi2=M_PI*M_PI;
//   dboot_t rl=MLep/MPS,rl2=rl*rl;
//   dboot_t bIR=1/(8*pi2)*(1+rl2)/(1-rl2)*log(rl2);
//   dboot_t b0=1/(16*pi2)*(-log(MW/MPS)+37/12.0+log(rl2)*(2*(1-3*rl2)+(1+rl2)*log(rl2))/(2*(1+rl2))+(1+rl2)/(1-rl2)*log(rl2)*(gammaeul-log(4*M_PI))+2*(1+rl2)*(k31+k32));
// }

//! compute the energy of a twisted mass quark
template <class T1,class T2>
T1 tm_quark_energy(T1 pi,T2 mass)
{
  T2 m2=mass*mass;
  T1 sinph=sin(pi/2);
  T1 sinph2=sinph*sinph;
  T1 sinph4=sinph2*sinph2;
  T1 p2=12*sinph2;
  T1 p4=12*sinph4;
  T1 four_sinh2_Eh=(m2+p2+p2*p2/4-p4)/(1+p2/2);
  
  return 2*asinh((T1)sqrt(four_sinh2_Eh/4));
}

//! compute the energy of a naive massless fermion
template <class T>
T naive_massless_quark_energy(const T &pi)
{
  T sinh2E=3*sqr(sin(pi));
  return asinh(sqrt(sinh2E));
}

template <class T>
T ratio_Wreg_leptonic_traces(const T &pi,const T &mass)
{
  T sinp=3*sqr(sin(pi));
  T sinhE=sinh(tm_quark_energy(pi,mass));
  
  return mass/(sinhE-sqrt(sinp));
}

//! compute the non-offshellness
double offshellness(double pi,double lep_mass,double mes_mass)
{
  double lep_energy=tm_quark_energy(pi,lep_mass);
  double neu_energy=naive_massless_quark_energy(pi);
  double err=lep_energy+neu_energy-mes_mass;
  
  return err;
}

//! find pi to put for a given process
double find_pi(double lep_mass,double mes_mass)
{return Brent_solve(bind(offshellness,_1,lep_mass,mes_mass),0,1);}

//! wrapper
template <class TS>
TS find_pi(double lep_mass,TS mes_mass)
{
  TS pi;
  for(size_t iel=0;iel<pi.size();iel++) pi[iel]=find_pi(lep_mass,mes_mass[iel]);
  return pi;
}

//! load or compute and write the z to compute FSE
djvec_t prepare_z_for_FSE(const string &path,const djack_t &jbetal)
{
  djvec_t jz(nZ_FSE);
  
  bool to_compute=true;
  
  //check that marker agrees
  if(file_exists(path))
    {
      cout<<"Z file "<<path<<" found, parsing"<<endl;
      raw_file_t fin(path,"r");
      djack_t jbetal_temp=fin.bin_read<djack_t>();
      to_compute=not equal(jbetal_temp.begin(),jbetal_temp.end(),jbetal.begin());
      cout<<"checking jbetal, new="<<smart_print(jbetal.ave_err())<<", stored="<<smart_print(jbetal_temp.ave_err())<<endl;
      if(not to_compute)
	{
	  cout<<"jbetal stored agree, reading z"<<endl;
	  fin.bin_read(jz);
	}
      else cout<<"jbetal stored does not agree"<<endl;
    }
  
  //if to compute, write it
  if(to_compute)
    {
      cout<<"Computing z, storing in "<<path<<endl;
      jz=zeta_FSE(jbetal);
      raw_file_t fout(path,"w");
      fout.bin_write(jbetal);
      fout.bin_write(jz);
    }
  
  return jz;
}

//! implements eq.50 of Silvano note BUT not dA/A
dboot_t Wreg1_contr(const dboot_t &a)
{
  const double uss=1/(16*sqr(M_PI));

  //tlSym photon
  // dboot_t Z1=uss*(5.0*log(a*MW)-5.056);
  // double Z2=uss*0.323;
  
  dboot_t Zmarci=uss*(-1.5-2.0*log(a*MW)-11.852);
  
  //pure Wilson photon
  dboot_t Z11=uss*(4.0*log(a*MW)-15.539)-0.5*Zmarci;
  double Z12=uss*0.536;
  
  // double Z21=uss*0.536;
  // dboot_t Z22=uss*(2.0*log(a*MW)-14.850)-0.5*Zmarci;
  
  //cout<<"Z1 (e2 included): "<<dboot_t(Z1*e2).ave_err()<<", Z2 (idem): "<<Z2*e2<<endl;
  
  //return 0.5*(Z11+Z12+Z21+Z22);
  return Z11+Z12;
}

//! (mixing and rotation op.2)
dboot_t Wreg2_contr(const dboot_t &a)
{
  const double uss=1/(16*sqr(M_PI));
  
  dboot_t Zmarci=uss*(-1.5-2.0*log(a*MW)-11.852);
  
  //pure Wilson photon
  dboot_t Z11=uss*(4.0*log(a*MW)-15.539)-0.5*Zmarci;
  dboot_t Z22=uss*(2.0*log(a*MW)-14.850)-0.5*Zmarci;
  
  return Z11-Z22;
}

//////////////////////////////////////////////////////////////////// cont chir extrap for hl //////////////////////////////////////////////////////


typedef vector<pair<size_t,double>> procs_t;

//! lepton energy
template <class Tpars>
Tpars elep(const double &MLep,const Tpars &M2PS)
{
  Tpars rl2=sqr(MLep)/M2PS;
  Tpars en=0.5*sqrt(M2PS)*(1.0+rl2);
  
  //cout<<"MLep: "<<MLep<<", M2PS: "<<M2PS<<", en: "<<en<<endl;
  
  return en;
}

//! finite size effects
// template <class Tpars,class Txi,class TL>
// Tpars FSE_corr_xi_hl(const Tpars &L2dep,const Tpars &L3dep,const double &MLep,const Txi &M2PS,const TL &L)
// {
//   return
//     L2dep/(M2PS*sqr(L))+
//     L3dep/sqr(Tpars(elep(MLep,M2PS)*L));
// }

//! finite size effects
template <class Tpars,class TL,class Tm>
Tpars FSE_corr_hl(const Tpars &L2dep,const Tpars &L3dep,const double &MLep,const Tpars &B0,const vector<Tm> &m,const procs_t &procs,const TL &L)
{
  Tpars out;
  out=0.0;
  
  for(auto &proc : procs)
    {
      size_t iproc=proc.first;
      size_t iQED_mes=iQED_mes_of_proc[iproc];
      size_t iq1=QED_mes_pars[iQED_mes].iq1;
      size_t iq2=QED_mes_pars[iQED_mes].iq2;
      Tpars M2PS=M2_fun(B0,m[iq1],m[iq2]);
      Tpars contr=proc.second*(L2dep/(M2PS*sqr(L))+L3dep/sqr(Tpars(elep(MLep,M2PS)*L)));
      out+=contr;
    }
  
  return out;
}

//! chiral behaviour
// template <class Tpars,class Txi>
// Tpars chir_corr_xi_hl(const Tpars &Kpi,const Tpars &K2pi,const Tpars &Z,const Txi &xi,const Txi &xi_s,const size_t iproc,const size_t chir_flag)
// {
//   Tpars out=Kpi*xi+K2pi*sqr(xi);
  
//   if(chir_flag==hl::chir::QUADRATICLOG)
//     {
//       switch(iproc)
// 	{
// 	case 0:
// 	  out+=(3.0-2.0*Z)*e2/(16*sqr(M_PI))*log(xi);
// 	  break;
// 	case 1:
// 	  out+=-(3.0-Z)*e2/(16*sqr(M_PI))*log(xi/xi_s);
// 	  break;
// 	default:
// 	  break;
// 	}
//     }
  
//   return out;
// }

//! chiral behaviour
template <class Tpars>
Tpars chir_corr_hl(const Tpars &Kpi,const Tpars &K2pi,const Tpars &Z,const Tpars &xi,const Tpars &xis,const procs_t &procs,const size_t chir_flag)
{
  const bool use_unquenched=true;
  
  Tpars out=Kpi*xi+K2pi*sqr(xi);
  //Tpars out=Kpi*(xis-xi)+K2pi*sqr(Tpars(xis-xi));
  
  if(chir_flag==hl::chir::QUADRATICLOG)
    for(auto &proc : procs)
      {
	size_t iproc=proc.first;
	Tpars chl;
	switch(iproc)
	  {
	  case 0:
	    if(use_unquenched) chl=(3.0-2.0*Z)*e2/(16*sqr(M_PI))*log(xi);
	    else               chl=(3.0-10.0*Z/9)*e2/(16*sqr(M_PI))*log(xi);
	    break;
	  case 1:
	    if(use_unquenched) chl=-Z*e2/(16*sqr(M_PI))*log(xi);
	    else               chl=-8.0*Z/9*e2/(16*sqr(M_PI))*log(xi);
	    break;
	  default:
	    chl=0;
	    break;
	  }
	out+=proc.second*chl;
      }
  
  return out;
}

//! ansatz fit
// template <class Tpars,class Txi,class Ta>
// Tpars cont_chir_ansatz_corr_xi_hl(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tpars &Z,const Txi &xi,const Txi &xi_s,const Txi &Mmes,const double MLep,const Ta &a,const Tpars &adep,double L,const Tpars &L2dep,const Tpars &L3dep,const procs_t &procs,const size_t chir_flag)
// {
//   Txi M2PS=sqr(Mmes);
//   Tpars res=C+
//     chir_corr_xi_hl(Kpi,K2pi,Z,xi,xi_s,procs,chir_flag)+
//     a*a*adep;
//   if(L>0) res+=FSE_corr_xi_hl(L2dep,L3dep,MLep,M2PS,L*a);
//   return res;
// }

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_corr_hl(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tpars &Z,const Tm &ml,const Tm &ms,const double MLep,const Ta &a,const Tpars &adep,double L,const Tpars &L2dep,const Tpars &L3dep,const procs_t &procs,const size_t chir_flag)
{
  Tpars xi=xi_fun(B0,ml,ml,f0);
  Tpars xi_s=xi_fun(B0,ml,ms,f0);
  return C+
    chir_corr_hl(Kpi,K2pi,Z,xi,xi_s,procs,chir_flag)+
    a*a*adep+
    FSE_corr_hl<Tpars,Ta,Tm>(L2dep,L3dep,MLep,B0,{ml,ms},procs,Ta(L*a));
}

template<class T>
void set_default_grace(const vector<T> &ext_data)
{
  //list of possible colors
  vector<grace::color_t> color_per_ib={grace::RED,grace::BLUE,grace::GREEN4};
  map<size_t,grace::symbol_t> symbol_per_L={{16,grace::STAR},{20,grace::CIRCLE},{24,grace::SQUARE},{32,grace::DIAMOND},{40,grace::TRIDOWN},{48,grace::TRIUP}};
  
  //make the list of volumes and beta
  set<pair<size_t,size_t>> list;
  for(size_t ib=0;ib<nbeta;ib++)
    {
      for(size_t idata=0;idata<ext_data.size();idata++)
	if(ext_data[idata].ib==ib)
	  list.insert(make_pair(ext_data[idata].ib,ext_data[idata].L));
    }
  
  //add the listed cols and symbols
  grace::default_symbol_scheme.clear();
  grace::default_color_scheme.clear();
  for(auto &p : list)
    {
      grace::default_color_scheme.push_back(color_per_ib[p.first]);
      grace::default_symbol_scheme.push_back(symbol_per_L[p.second]);
    }
  //add c.l
  grace::default_color_scheme.push_back(grace::VIOLET);
  grace::default_symbol_scheme.push_back(grace::TRILEFT);
}

//! perform the fit to the continuum limit of correction of process
// dboot_t cont_chir_fit_corr_xi_hl(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_xi_data_t> &ext_data,const dboot_t &xi_phys,const dboot_t &xi_s_phys,const double MMes_phys,const double &MLep,const string &path,const size_t iproc,const size_t isyst,const bool cov_flag,const vector<string> &beta_list)
// {
//   using namespace hl;
  
//   //set_printlevel(3);
  
//   boot_fit_t boot_fit;
//   size_t nbeta=a.size();
//   cont_chir_fit_xi_pars_t pars(nbeta);
  
//   //guesses
//   ave_err_t L2dep_guess;
//   ave_err_t L3dep_guess;
//   const ave_err_t C_guess[2]={{0.021,0.001},{-0.012,0.001}};
//   const ave_err_t KPi_guess(-0.34,0.01);
//   const ave_err_t K2Pi_guess(1.0,0.5);
//   const ave_err_t Z_guess={0.658,0.040};
//   const ave_err_t adep_guess={0.003,0.006};
//   const ave_err_t adep_ml_guess={0,0.001};
  
//   const size_t FSE_flag=FSE::variations[case_of<c_FSE>(isyst)];
//   switch(FSE_flag)
//     {
//       using namespace FSE;
//     case WITHSTDEP:
//       L2dep_guess=ave_err_t(-0.16,0.1);
//       L3dep_guess=ave_err_t(0.04,0.03);
//       break;
//     case NOSMALLVOL:
//       L2dep_guess=ave_err_t(-0.068,0.09);
//       L3dep_guess=ave_err_t(-0.03,0.03);
//       break;
//     case NOSTDEP:
//       L2dep_guess=ave_err_t(0.0,0.1);
//       L3dep_guess=ave_err_t(0.0,0.1);
//     }
  
//   //set parameters
//   pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
//   pars.iC=boot_fit.add_fit_par(pars.C,"C_guess",C_guess[iproc].ave(),C_guess[iproc].err());
//   pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
//   pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
//   pars.iL2dep=boot_fit.add_fit_par(pars.L2dep,"L2dep",L2dep_guess.ave(),L2dep_guess.err());
//   pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave(),L3dep_guess.err());
//   pars.iKK=boot_fit.add_self_fitted_point(pars.KK,"Z",Z_guess);
//   //boot_fit.fix_par_to(pars.iL3dep,0.0);
  
//   //set FSE pars
//   switch(FSE_flag)
//     {
//       using namespace FSE;
//     case NOSTDEP:
//       boot_fit.fix_par_to(pars.iL2dep,0.0);
//       boot_fit.fix_par_to(pars.iL3dep,0.0);
//       break;
//     case WITHSTDEP:
//     case NOSMALLVOL:
//       break;
//     }
  
//   //set cont limit pars
//   const size_t cont_flag=cont::variations[case_of<c_cont>(isyst)];
//   boot_fit.fix_par_to(pars.iadep_xi,0.0);
//   switch(cont_flag)
//     {
//       using namespace cont;
//     case(CONSTANT):
//       boot_fit.fix_par_to(pars.iadep,0.0);
//     break;
//     case(LINEAR):
//       break;
//     }
  
//   //set chir limit pars
//   const size_t chir_flag=chir::variations[case_of<c_chir>(isyst)];
//   switch(chir_flag)
//     {
//       using namespace chir;
//     case(QUADRATICLOG):
//       //boot_fit.fix_par_to(pars.iK2Pi,0.0);
//       break;
//     case(QUADRATIC):
//       //boot_fit.fix_par_to(pars.iKK,0.0);
//       break;
//     }
  
//   cont_chir_fit_xi_minimize(ext_data,pars,boot_fit,0.0,0.0,[iproc,chir_flag](const vector<double> &p,const cont_chir_fit_xi_pars_t &pars,double xi,double xi_s,double MMes,double MLep,double ac,double L)
// 			    {return cont_chir_ansatz_corr_xi_hl(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],p[pars.iKK],xi,xi_s,MMes,MLep,ac,p[pars.iadep],L,p[pars.iL2dep],p[pars.iL3dep],iproc,chir_flag);}
// 			 ,cov_flag);
  
//   dboot_t phys_res=cont_chir_ansatz_corr_xi_hl(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,xi_phys,xi_s_phys,dboot_t(MMes_phys),MLep,1,dboot_t(pars.adep*0.0),-1,(dboot_t)(0.0*pars.L2dep),dboot_t(0.0*pars.L3dep),iproc,chir_flag);
//   cout<<"result: "<<phys_res.ave_err()<<", "<<phys_res[0]<<endl;
  
//   //bool include_small_vol=(FSE::variations[case_of<c_FSE>(isyst)]!=hl::FSE::NOSTDEP);
//   //bool include_coarse=(cont::variations[case_of<c_cont>(isyst)]!=hl::cont::CONSTANT);
  
//   set_default_grace(ext_data);
  
//   const string yaxis_title="$$\\delta m_l^{ren}";
//   plot_chir_fit_xi(path,ext_data,pars,
// 		[&pars,&MLep,&xi_s_phys,&MMes_phys,iproc,chir_flag]
// 		(double x,size_t ib)
// 		{return cont_chir_ansatz_corr_xi_hl<double,double,double>
// 		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),pars.KK.ave(),x,xi_s_phys.ave(),MMes_phys,MLep,pars.fit_a[ib].ave(),pars.adep.ave(),
// 		     inf_vol,0.0*pars.L2dep.ave(),0.0*pars.L3dep.ave(),iproc,chir_flag);},
// 		   bind(cont_chir_ansatz_corr_xi_hl<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,_1,xi_s_phys.ave(),MMes_phys,MLep,1.0,dboot_t(pars.adep*0.0),
// 			-1,dboot_t(0.0*pars.L2dep),dboot_t(0.0*pars.L3dep),iproc,chir_flag),
// 		[&ext_data,&pars,&MLep]
// 		(size_t idata,bool without_with_fse,size_t ib)
// 		{
// 		  dboot_t a=pars.fit_a[ib];
// 		  dboot_t z=pars.fit_z[ib];
// 		  dboot_t xi=ext_data[idata].xi;
// 		  dboot_t xi_s=ext_data[idata].xi_s;
// 		  dboot_t MMes=ext_data[idata].MMes;
// 		  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_corr_xi_hl(pars.L2dep,pars.L3dep,MLep,sqr(MMes),ext_data[idata].L*a));},
// 		   xi_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
//   return phys_res;
// }

//! perform the fit to the continuum limit of correction of process
dboot_t cont_chir_fit_corr_hl(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const double &MLep,const string &path,const procs_t &procs,const size_t isyst,const bool cov_flag,const vector<string> &beta_list,const size_t istudy,const vector<hl::FSE::SUB_STDEP::t> FSE_sub_variations,const vector<hl::FSE::FIT_STDEP::t> &FSE_fit_variations)
{
  using namespace hl;
  
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  //guesses
  ave_err_t L2dep_guess;
  ave_err_t L3dep_guess;
  ave_err_t C_guess;
  ave_err_t KPi_guess;
  ave_err_t K2Pi_guess;
  const ave_err_t Z_guess={0.658,0.040};
  ave_err_t adep_guess;
  const ave_err_t adep_ml_guess={0.0,0.001};
  
  const size_t FSE_FIT_flag=FSE_fit_variations[case_of<c_FSE>(isyst)];
  const size_t FSE_SUB_flag=FSE_sub_variations[case_of<c_FSE>(isyst)];
  const size_t chir_flag=chir::variations[case_of<c_chir>(isyst)];
  
  switch(istudy)
    {
    case STUDY_PI:
      C_guess={0.021,0.001};
      KPi_guess={-0.18,0.2};
      K2Pi_guess={1.0,0.5};
      adep_guess={-0.0003,0.006};
      
      //we always fit FSE
      L2dep_guess=ave_err_t(0.16,0.2);
      L3dep_guess=ave_err_t(-0.04,0.03);
      break;
      
    case STUDY_K_M_PI:
      adep_guess={-0.005,0.003};
      
      switch(chir_flag)
	{
	  using namespace chir;
	case QUADRATIC:
	  C_guess={-0.014,0.002};
	  KPi_guess={0.22,0.06};
	  K2Pi_guess={-1.3,0.5};
	  break;
	case QUADRATICLOG:
	  C_guess={-0.023,0.002};
	  KPi_guess={0.30,0.08};
	  K2Pi_guess={-1.7,0.7};
	  break;
	}
      
      switch(FSE_FIT_flag)
	{
	  using namespace FSE::FIT_STDEP;
	case YES:
	  L2dep_guess=ave_err_t(-0.27,0.1);
	  L3dep_guess=ave_err_t(0.1,0.03);
	  break;
	case NO:
	  L2dep_guess=ave_err_t(0.0,0.1);
	  L3dep_guess=ave_err_t(0.0,0.1);
	  break;
	}
      break;
      
    case STUDY_K:
      C_guess={0.003,0.001};
      KPi_guess={0.00,0.03};
      K2Pi_guess={0.0,0.1};
      adep_guess={0.003,0.006};
      
      //we always fit FSE
      //Z_fact=1.0
      switch(FSE_SUB_flag)
	{
	  using namespace FSE::SUB_STDEP;
	case YES:
	  L2dep_guess=ave_err_t(1.1,8.7);
	  L3dep_guess=ave_err_t(-0.4,2.3);
	  break;
	case NO:
	  L2dep_guess=ave_err_t(0.5,2.0);
	  L3dep_guess=ave_err_t(-0.2,0.5);
	  break;
	}
      //Z_fact=0.75
      // switch(FSE_SUB_flag)
      // 	{
      // 	  using namespace FSE::SUB_STDEP;
      // 	case YES:
      // 	  L2dep_guess=ave_err_t(-1.5,2.0);
      // 	  L3dep_guess=ave_err_t(-0.04,0.03);
      // 	  break;
      // 	case NO:
      // 	  L2dep_guess=ave_err_t(0.16,0.1);
      // 	  L3dep_guess=ave_err_t(-0.04,0.03);
      // 	  break;
      // 	}
      break;
    }
  
  //set parameters
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C_guess",C_guess.ave(),C_guess.err());
  cout<<"C_guess "<<C_guess.ave()<<" "<<C_guess.err()<<endl;
  cout<<"L2dep_guess "<<L2dep_guess.ave()<<" "<<L2dep_guess.err()<<endl;
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave(),KPi_guess.err());
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",K2Pi_guess.ave(),K2Pi_guess.err());
  pars.iL2dep=boot_fit.add_fit_par(pars.L2dep,"K2",L2dep_guess.ave(),L2dep_guess.err());
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"K2ell",L3dep_guess.ave(),L3dep_guess.err());
  pars.iKK=boot_fit.add_self_fitted_point(pars.KK,"Z",Z_guess);
  //boot_fit.fix_par_to(pars.iL3dep,0.0);
  
  //set FSE pars
  switch(FSE_FIT_flag)
    {
      using namespace FSE::FIT_STDEP;
    case NO:
      cout<<"Fixing Ldep"<<endl;
      boot_fit.fix_par_to(pars.iL2dep,0.0);
      boot_fit.fix_par_to(pars.iL3dep,0.0);
      break;
    case YES:
      break;
    }
  
  //set cont limit pars
  const size_t cont_flag=cont::variations[case_of<c_cont>(isyst)];
  boot_fit.fix_par_to(pars.iadep_ml,0.0);
  switch(cont_flag)
    {
      using namespace cont;
    case(CONSTANT):
      boot_fit.fix_par_to(pars.iadep,0.0);
    break;
    case(LINEAR):
      break;
    }
  
  //set chir limit pars
  switch(chir_flag)
    {
      using namespace chir;
    case(QUADRATICLOG):
      //boot_fit.fix_par_to(pars.iK2Pi,0.0);
      break;
    case(QUADRATIC):
      //boot_fit.fix_par_to(pars.iKK,0.0);
      break;
    }
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[procs,chir_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double MLep,double ac,double L)
                         {return cont_chir_ansatz_corr_hl(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],p[pars.iKK],ml,ms,MLep,ac,p[pars.iadep],L,p[pars.iL2dep],p[pars.iL3dep],procs,chir_flag);}
                         ,cov_flag);
  
  double a_cont=1.0e-5;
  dboot_t phys_res=cont_chir_ansatz_corr_hl(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,ml_phys,ms_phys,MLep,a_cont,pars.adep,inf_vol,(dboot_t)(0.0*pars.L2dep),(dboot_t)(0.0*pars.L3dep),procs,chir_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  //bool include_small_vol=(FSE::variations[case_of<c_FSE>(isyst)]!=hl::FSE::NOSTDEP);
  //bool include_coarse=(cont::variations[case_of<c_cont>(isyst)]!=hl::cont::CONSTANT);
  
  set_default_grace(ext_data);
  
  const string yaxis_title="$$\\delta m_l^{ren}";
  plot_chir_fit(path,ext_data,pars,
                [&pars,&MLep,&ms_phys,procs,chir_flag]
                (double x,size_t ib)
                {return cont_chir_ansatz_corr_hl<double,double,double>
                    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),pars.KK.ave(),x,ms_phys.ave(),MLep,pars.fit_a[ib].ave(),pars.adep.ave(),
                     inf_vol,0.0*pars.L2dep.ave(),0.0*pars.L3dep.ave(),procs,chir_flag);},
                bind(cont_chir_ansatz_corr_hl<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,_1,ms_phys.ave(),MLep,a_cont,pars.adep,
                     inf_vol,(dboot_t)(0.0*pars.L2dep),(dboot_t)(0.0*pars.L3dep),procs,chir_flag),
                [&ext_data,&pars,&MLep,&procs]
                (size_t idata,bool without_with_fse,size_t ib)
                {
                  dboot_t a=pars.fit_a[ib];
                  dboot_t z=pars.fit_z[ib];
                  dboot_t ml=ext_data[idata].aml/a/z;
                  dboot_t ms=ext_data[idata].ams/a/z;
                  return dboot_t(ext_data[idata].wfse-without_with_fse*FSE_corr_hl<dboot_t,dboot_t,dboot_t>(pars.L2dep,pars.L3dep,MLep,pars.fit_B0,{ml,ms},procs,ext_data[idata].L*a));},
                ml_phys,phys_res,yaxis_title,beta_list,ind_syst.descr(isyst));
  
  //! output file for A40 data
  vector<size_t> syst_comp=ind_syst(isyst);
  size_t iFSE=syst_comp[c_FSE];
  if(accumulate(syst_comp.begin(),syst_comp.end(),-iFSE)==0)
    {
      grace_file_t A40_XX_file_FSE_full_sub;
      size_t iFSE=case_of<c_FSE>(isyst);
      A40_XX_file_FSE_full_sub.open(combine("plots_hl/A40_study%zu_FSEord%zu_fullsub.xmg",istudy,iFSE));
      
      dboot_t ml;
      
      //! two sets of data, with and without FSE
      for(size_t without_with=0;without_with<=1;without_with++)
	{
	  A40_XX_file_FSE_full_sub.new_data_set();
	  for(size_t iens=0;iens<nens_used;iens++)
	    {
	      const ens_pars_t &ens=ens_pars[iens];
	      const size_t ib=ens.ib;
	      if(fabs(ens_pars[iens].aml-0.0040)<1e-6)
		{
		  const dboot_t &a=pars.fit_a[ib];
		  const dboot_t &z=pars.fit_z[ib];
		  ml=ext_data[iens].aml/a/z;
		  dboot_t ms=ext_data[iens].ams/a/z;
		  const dboot_t &out=ext_data[iens].wfse-without_with*
		    FSE_corr_hl(pars.L2dep,pars.L3dep,MLep,pars.fit_B0,vector<dboot_t>{ml,ms},procs,ext_data[iens].L*a);
		  A40_XX_file_FSE_full_sub.write_ave_err(sqr(1.0/ens.L),out.ave_err());
		}
	    }
	}
      
      //fit band
      const double x_min=sqr(0.001),x_max=5e-3;
      A40_XX_file_FSE_full_sub.write_polygon(
		    [&pars,&MLep,&ms_phys,&ml,procs,chir_flag]
		    (double x)
		    {return cont_chir_ansatz_corr_hl<dboot_t,dboot_t,dboot_t>
			(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,ml,ms_phys,MLep,pars.fit_a[0],pars.adep,
			 1/sqrt(x),pars.L2dep,pars.L3dep,procs,chir_flag);},
		    x_min,x_max);
      
      //infinite volume extrapolate
      dboot_t out=cont_chir_ansatz_corr_hl(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.KK,ml,ms_phys,MLep,pars.fit_a[0],pars.adep,inf_vol,(dboot_t)(0.0*pars.L2dep),(dboot_t)(0.0*pars.L3dep),procs,chir_flag);
      A40_XX_file_FSE_full_sub.write_constant_band(x_min,x_max,out);
    }
  
  return phys_res;
}

//! compute the correction to the process
dbvec_t compute_corr(size_t iproc,const int &include_stong_IB)
{
  //! data to extrapolate
  dbvec_t tot_corr_all(hl::ind_an_ens_FSEmax_frange.max());
  
  //open a table for FSE contribution
  string FSE_tab_path="tables/FSE_proc"+to_string(iproc)+".txt";
  ofstream FSE_tab(FSE_tab_path);
  if(not FSE_tab.good()) CRASH("unable to open %s",FSE_tab_path.c_str());
  
  //open a table for Wreg1 contribution
  string Wreg1_contr_tab_path="tables/Wreg1_contr_proc"+to_string(iproc)+".txt";
  ofstream Wreg1_contr_tab(Wreg1_contr_tab_path);
  if(not Wreg1_contr_tab.good()) CRASH("unable to open %s",Wreg1_contr_tab_path.c_str());

  //open a table for Wreg2 contribution
  string Wreg2_contr_tab_path="tables/Wreg2_contr_proc"+to_string(iproc)+".txt";
  ofstream Wreg2_contr_tab(Wreg2_contr_tab_path);
  if(not Wreg2_contr_tab.good()) CRASH("unable to open %s",Wreg2_contr_tab_path.c_str());
  
  //open a table for QED contribution
  string qed_corr_tab_path="tables/QED_corr_proc"+to_string(iproc)+".txt";
  ofstream qed_corr_tab(qed_corr_tab_path);
  if(not qed_corr_tab.good()) CRASH("unable to open %s",qed_corr_tab_path.c_str());
  
  //open a table for bc
  string bc_tab_path="tables/bc_table"+to_string(iproc)+".txt";
  static ofstream bc_tab(bc_tab_path);
  if(not bc_tab.good()) CRASH("unable to open %s",bc_tab_path.c_str());
  
  const size_t ilep=iMLep_of_proc[iproc];
  const size_t iQED_mes=iQED_mes_of_proc[iproc];
  const size_t iQCD_mes=QED_mes_pars[iQED_mes].iQCD;
  
  //! output file for A40 data
  grace_file_t A40_XX_file(combine("plots_hl/A40_proc%zu.xmg",iproc));
  vector<grace_file_t> A40_XX_file_FSE_sub(hl::FSE_max_orders.size());
  for(size_t iFSE_max=0;iFSE_max<hl::FSE_max_orders.size();iFSE_max++)
    A40_XX_file_FSE_sub[iFSE_max].open(combine("plots_hl/A40_proc%zu_FSEord%zu.xmg",iproc,hl::FSE_max_orders[iFSE_max]));
  
  //! loop over analysis
  for(size_t iens=0;iens<nens_used;iens++)
    {
      const size_t ind_QCD=ind_ens_QCD_mes({iens,iQCD_mes});
      const size_t ind_QED=ind_ens_QED_mes({iens,iQED_mes});
      const ens_pars_t &ens=ens_pars[iens];
      const double aMLep=ens.aMLep[ilep];
      const djack_t jpi=find_pi(aMLep,jaM[ind_QCD]);
      
      const djack_t jbetal=sqrt(3)*jpi/tm_quark_energy(jpi,aMLep);
      djvec_t jz(nZ_FSE);
      if(iproc<2) jz=prepare_z_for_FSE(ens.path+"/z_FSE_proc"+to_string(iproc)+".dat",jbetal);
      else        jz=0.0;
      
      const size_t ib=ens.ib;
      const size_t ind_proc=ind_ens_proc({iens,iproc});
      
      //test
      const djvec_t loop=jLO_A_bare[ind_proc]/jAP_LO_exp_removed[ind_QCD];
      jAP_LO_exp_removed[ind_QCD].ave_err().write(combine("%s/plots_hl/loop_den_iproc%zu.xmg",ens.path.c_str(),iproc));
      loop.ave_err().write(combine("%s/plots_hl/loop_iproc%zu.xmg",ens.path.c_str(),iproc));
      
      for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
	{
	  const size_t ind_an_ens=ind_adml({input_an_id,iens});
	  
	  //compute the nasty diagram
	  const dbvec_t LO=Zv[ib]*dbvec_t(bi,-jLO_A_bare[ind_proc]);
	  const dbvec_t QED=Zv[ib]*dbvec_t(bi,-jQED_A_bare[ind_proc])+dbvec_t(bi,jQED_V_bare[ind_proc])*Za[ib]; //minus because V-A
	  LO.ave_err().write(combine("%s/plots_hl/LO_iproc%zu.xmg",ens.path.c_str(),iproc));
	  QED.ave_err().write(combine("%s/plots_hl/QED_iproc%zu.xmg",ens.path.c_str(),iproc));
	  
	  //check bc put in the simulation
	  const double pi_bare=find_pi(aMLep,ens.aMMes[iQCD_mes]);
	  const double bc=pi_bare*ens.L/M_PI;
	  bc_tab<<"bc put in the simulation for ensemble "<<ens.path<<": "<<bc<<endl;
	  
	  //cout<<"betal for ensemble "<<ens.path<<": "<<betal.ave_err()<<endl;
	  const dboot_t a=1/lat_par[input_an_id].ainv[ib];
	  const dboot_t L=ens.L*a;
	  const dboot_t MLep=aMLep/a;
	  const dboot_t Mmes=dboot_t(bi,jaM[ind_QCD])/a;
	  
	  //! FSE contribution: this contains already a 2
	  const dbvec_t FSE_contr=FSE_corr(MLep,Mmes,z0,dbvec_t(bi,jz),L,hl::FSE_max_orders)*e2;
	  FSE_tab<<"Ensemble: "<<ens.path<<", proc: "<<iproc<<", input_an_id: "<<input_an_id<<endl;
	  for(size_t i=0;i<FSE_contr.size();i++) FSE_tab<<" "<<i<<" "<<FSE_contr[i].ave_err()<<endl;
	  
	  //! contribution due to W reg1 (2*e2 added)
	  const dboot_t W1_contr=2*Wreg1_contr(a)*e2;
	  Wreg1_contr_tab<<"Ensemble: "<<ens.path<<", proc: "<<iproc<<", input_an_id: "<<input_an_id<<", Wreg1_contr to amplitude: "<<
	    smart_print(W1_contr.ave_err())<<endl;

	  //! contribution due to W reg2 (2*e2 added)
	  const dboot_t W2_contr=2*Wreg2_contr(a)*e2*(Za[ib]-Zv[ib])/(2*Zv[ib]);
	  Wreg2_contr_tab<<"Ensemble: "<<ens.path<<", proc: "<<iproc<<", input_an_id: "<<input_an_id<<", Wreg2_contr to amplitude: "<<
	    smart_print(W2_contr.ave_err())<<endl;
	  
	  //! maximum energy for emitted photon
	  const dboot_t DeltaE=Mmes*(1-sqr((dboot_t)(MLep[ilep]/Mmes)))/2.0;
	  cout<<"DeltaE: "<<DeltaE.ave_err()<<endl;
	  
	  //! extract dA/A from nasty diagram ratio
	  dbvec_t rat_ext=QED/LO;
	  rat_ext[rat_ext.size()-1]=rat_ext[0]=0.0; //set to zero the contact term
	  const size_t itint=QCD_mes_pars[iQCD_mes].itint;
	  for(size_t ifrange=0;ifrange<nfit_range_variations;ifrange++)
	    {
	      const size_t rvar=frange_var[ifrange],tin=std::min(ens.tmin[itint],ens.T/4-1-rvar);
	      const size_t tmin=tin+rvar,tmax=ens.T/2-tin-rvar;
	      const string dA_fr_A_path=combine("%s/plots_hl/QED_LO_ratio_iproc%zu_ian%zu_frange%zu.xmg",ens.path.c_str(),iproc,input_an_id,ifrange);
	      dboot_t dA_fr_A=constant_fit(rat_ext,tmin,tmax,dA_fr_A_path);
	      
	      //compute the internal+external contribution
	      const dboot_t external=2.0*dA_fr_A*e2;
	      const dboot_t dZA_QED_rel=dboot_t(bi,jDZA_QED_rel[ind_QED])*e2;
	      const dboot_t dZA_MASS_rel=dboot_t(bi,jDZA_MASS_rel[ind_QED])*adml_bare[ind_an_ens];
	      const dboot_t dM_QED_rel=dboot_t(bi,jDM_QED[ind_QED]/jaM[ind_QCD])*e2;
	      const dboot_t dM_MASS_rel=dboot_t(bi,jDM_MASS[ind_QED]/jaM[ind_QCD])*adml_bare[ind_an_ens];
	      const dboot_t internal_QED=2.0*dZA_QED_rel;
	      const dboot_t internal_MASS=2.0*dZA_MASS_rel*include_stong_IB;
	      const dboot_t rate_QED_mass=-2.0*dM_QED_rel; //to be SUBTRACTED
	      const dboot_t rate_MASS_mass=-2.0*dM_MASS_rel*include_stong_IB; //and this as well
	      const double marc_sirl=e2/(2*sqr(M_PI))*log(MZ/MW);
	      const dboot_t rate_pt=Gamma_pt(MLep,Mmes,DeltaE)*e2; //only e2

	      const double Z_fact=1.0;
	      
	      const dboot_t tot_but_FSE=external+internal_QED+internal_MASS+(W1_contr+W2_contr)*Z_fact+rate_QED_mass+rate_MASS_mass+rate_pt+marc_sirl;
	      if(ifrange==0 and input_an_id==0 and ens.ib==0 and fabs(ens.aml-0.0040)<1e-6)
		A40_XX_file.write_ave_err(ens.L,tot_but_FSE.ave_err());
	      
	      for(size_t iFSE_max=0;iFSE_max<hl::FSE_max_orders.size();iFSE_max++)
		{
		  const dboot_t tot_corr=tot_but_FSE-FSE_contr[iFSE_max];
		  if(ifrange==0 and input_an_id==0 and ens.ib==0 and fabs(ens.aml-0.0040)<1e-6)
		    A40_XX_file_FSE_sub[iFSE_max].write_ave_err(ens.L,tot_corr.ave_err());
		  
		  qed_corr_tab<<"Ensemble: "<<ens.path<<", proc: "<<iproc<<", input_an_id: "<<input_an_id<<", ifrange: "<<ifrange<<", iFSE_max: "<<iFSE_max<<
		    ",\n external: "<<smart_print(external.ave_err())<<
		    ",\n internal QED: "<<smart_print(internal_QED.ave_err())<<
		    ",\n internal MASS: "<<smart_print(internal_MASS.ave_err())<<
		    ",\n W1 contr: "<<smart_print(W1_contr.ave_err())<<
		    ",\n W2 contr: "<<smart_print(W2_contr.ave_err())<<
		    ",\n FSE contr: "<<smart_print(FSE_contr[iFSE_max].ave_err())<<
		    ",\n rate QED mass: "<<smart_print(rate_QED_mass.ave_err())<<
		    ",\n rate MASS mass: "<<smart_print(rate_MASS_mass.ave_err())<<
		    ",\n  (mass): "<<smart_print(jaM[ind_QCD].ave_err())<<
		    ",\n rate pt: "<<smart_print(rate_pt.ave_err())<<
		    ",\n marc sirl: "<<marc_sirl<<
		    ",\n tot corr to rate: "<<smart_print(tot_corr.ave_err())<<endl<<endl;
		  
		  cout<<"Tot: "<<tot_corr.ave_err()<<endl;
		  size_t i=hl::ind_an_ens_FSEmax_frange({input_an_id,iens,iFSE_max,ifrange});
		  tot_corr_all[i]=tot_corr;
		}
	    }
	}
    }
  
  return tot_corr_all;
}

//! extrapolate a single corr
void extrapolate_corr(const dbvec_t &tot_corr_all,const procs_t &procs,const size_t istudy,const size_t ilep,
		      const vector<hl::FSE::SUB_STDEP::t> FSE_sub_variations,const vector<hl::FSE::FIT_STDEP::t> &FSE_fit_variations)
{
  dbvec_t res(hl::ind_syst.max());
  
  //! perform fit
  for(size_t isyst=0;isyst<hl::ind_syst.max();isyst++)
    {
      cout<<"------------------------------------------------- "<<isyst<<" -------------------------------------------------"<<endl;
      
      cout<<hl::ind_syst.descr(isyst)<<endl;
      const size_t iFSE_variation=hl::case_of<hl::c_FSE>(isyst);
      const hl::FSE::FIT_STDEP::t FSE_fit=FSE_fit_variations[iFSE_variation];
      const hl::FSE::SUB_STDEP::t FSE_sub=FSE_sub_variations[iFSE_variation];
      cout<<"CONT: "<<hl::cont::tag[hl::case_of<hl::c_cont>(isyst)]<<endl;
      cout<<"CHIR: "<<hl::chir::tag[hl::case_of<hl::c_chir>(isyst)]<<endl;
      cout<<"FSE_FIT: "<<hl::FSE::FIT_STDEP::tag[FSE_fit]<<endl;
      cout<<"FSE_SUB: "<<hl::FSE::SUB_STDEP::tag[FSE_sub]<<endl;
      
      const size_t input_an_id=hl::case_of<hl::c_input>(isyst);
      const size_t ifrange=hl::case_of<hl::c_frange>(isyst);
      const size_t chir_FLAG=hl::chir::variations[hl::case_of<hl::c_chir>(isyst)];
      //const size_t cont_FLAG=hl::cont::variations[hl::case_of<hl::c_cont>(isyst)];
      const bool use_cov=false;
      
      //! add to the fit
#ifdef XI
      vector<cont_chir_fit_xi_data_t> fit_data;
#else
      vector<cont_chir_fit_data_t> fit_data;
#endif
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  const ens_pars_t &ens=ens_pars[iens];
	  const size_t idata=hl::ind_an_ens_FSEmax_frange({input_an_id,iens,FSE_sub,ifrange});
	  
	  //check if to include
	  bool include=true;
	  //include&=ens.use_for_L;
	  if(chir_FLAG==hl::chir::QUADRATICLOG) include&=ens.use_for_chir;
	  //if(cont_FLAG==hl::cont::CONSTANT) include&=ens.use_for_a;
	  
	  if(include)
	    {
	      dboot_t aMaux;
	      aMaux=ens.aMLep[ilep];
#ifdef XI
	      const size_t iens_QCD_mes=ind_ens_QCD_mes({iens,iQCD_mes});
	      dboot_t xi=dboot_t(bi,jxi[iens_QCD_mes]);
	      dboot_t xi_s=dboot_t(bi,jxi[iens_QCD_mes]);
	      dboot_t MMes=dboot_t(bi,jaM[iens_QCD_mes]);
	      fit_data.push_back(cont_chir_fit_xi_data_t(xi,xi_s,MMes,aMaux,ens.ib,ens.L,tot_corr_all[idata],tot_corr_all[idata]));
#else
              fit_data.push_back(cont_chir_fit_data_t(ens.aml,ens.ams,aMaux,ens.ib,ens.L,tot_corr_all[idata],tot_corr_all[idata]));
#endif
	    }
	}
      
      const string cc_path=combine("plots_hl/cont_chir_study%zu_isyst%zu.xmg",istudy,isyst);
      
#ifdef XI
      const dboot_t xi_phys=sqr(dboot_t(MP0/(4*M_PI*fP0)));
      const dboot_t xi_s_phys=sqr(dboot_t(MK0/(4*M_PI*fK0)));
      const double MMes_phys((iproc==0)?MP0:MK0);
      res[isyst]=cont_chir_fit_corr_xi_hl(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,fit_data,xi_phys,xi_s_phys,MMes_phys,MLep[iproc],cc_path,iproc,isyst,use_cov,beta_list);
#else
      res[isyst]=cont_chir_fit_corr_hl(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,fit_data,lat_par[input_an_id].ml,lat_par[input_an_id].ms,MLep[ilep],cc_path,procs,isyst,use_cov,beta_list,istudy,FSE_sub_variations,FSE_fit_variations);
#endif
    }
  
  perform_analysis(res,hl::ind_syst,"Res");
}

void test_factorization(size_t iproc)
{
  set_default_grace(ens_pars);

  load_all_hl({1,0},{1,0});
  const vector<djvec_t> jLO_A_bare_00=jLO_A_bare;
  const vector<djvec_t> jQED_qA0_lV0_bare_00=jQED_qA0_lV0_bare;
  const vector<djvec_t> jQED_qV0_lV0_bare_00=jQED_qV0_lV0_bare;
  const vector<djvec_t> jQED_qAi_lVi_bare_00=jQED_qAi_lVi_bare;
  const vector<djvec_t> jQED_qVi_lVi_bare_00=jQED_qVi_lVi_bare;

  load_all_hl({1,0},{0,1});
  const vector<djvec_t> jLO_A_bare_01=jLO_A_bare;
  const vector<djvec_t> jQED_qA0_lV0_bare_01=jQED_qA0_lV0_bare;
  const vector<djvec_t> jQED_qV0_lV0_bare_01=jQED_qV0_lV0_bare;
  const vector<djvec_t> jQED_qAi_lVi_bare_01=jQED_qAi_lVi_bare;
  const vector<djvec_t> jQED_qVi_lVi_bare_01=jQED_qVi_lVi_bare;

  load_all_hl({0,1},{1,0});
  const vector<djvec_t> jLO_A_bare_10=jLO_A_bare;
  const vector<djvec_t> jQED_qA0_lV0_bare_10=jQED_qA0_lV0_bare;
  const vector<djvec_t> jQED_qV0_lV0_bare_10=jQED_qV0_lV0_bare;
  const vector<djvec_t> jQED_qAi_lVi_bare_10=jQED_qAi_lVi_bare;
  const vector<djvec_t> jQED_qVi_lVi_bare_10=jQED_qVi_lVi_bare;

  load_all_hl({0,1},{0,1});
  const vector<djvec_t> jLO_A_bare_11=jLO_A_bare;
  const vector<djvec_t> jQED_qA0_lV0_bare_11=jQED_qA0_lV0_bare;
  const vector<djvec_t> jQED_qV0_lV0_bare_11=jQED_qV0_lV0_bare;
  const vector<djvec_t> jQED_qAi_lVi_bare_11=jQED_qAi_lVi_bare;
  const vector<djvec_t> jQED_qVi_lVi_bare_11=jQED_qVi_lVi_bare;

  const vector<djvec_t> jLO_A_bare_1=(jLO_A_bare_00+jLO_A_bare_11)/2.0;
  const vector<djvec_t> jQED_qA0_lV0_bare_1=(jQED_qA0_lV0_bare_00+jQED_qA0_lV0_bare_11)/2.0;
  const vector<djvec_t> jQED_qV0_lV0_bare_1=(jQED_qV0_lV0_bare_00+jQED_qV0_lV0_bare_11)/2.0;
  const vector<djvec_t> jQED_qAi_lVi_bare_1=(jQED_qAi_lVi_bare_00+jQED_qAi_lVi_bare_11)/2.0;
  const vector<djvec_t> jQED_qVi_lVi_bare_1=(jQED_qVi_lVi_bare_00+jQED_qVi_lVi_bare_11)/2.0;
  const vector<djvec_t> jQED_A_bare_1=jQED_qA0_lV0_bare_1+jQED_qAi_lVi_bare_1;
  const vector<djvec_t> jQED_V_bare_1=jQED_qV0_lV0_bare_1+jQED_qVi_lVi_bare_1;

  const vector<djvec_t> jLO_A_bare_0=(jLO_A_bare_01+jLO_A_bare_10)/2.0;
  const vector<djvec_t> jQED_qA0_lV0_bare_0=(jQED_qA0_lV0_bare_01+jQED_qA0_lV0_bare_10)/2.0;
  const vector<djvec_t> jQED_qV0_lV0_bare_0=(jQED_qV0_lV0_bare_01+jQED_qV0_lV0_bare_10)/2.0;
  const vector<djvec_t> jQED_qAi_lVi_bare_0=(jQED_qAi_lVi_bare_01+jQED_qAi_lVi_bare_10)/2.0;
  const vector<djvec_t> jQED_qVi_lVi_bare_0=(jQED_qVi_lVi_bare_01+jQED_qVi_lVi_bare_10)/2.0;
  const vector<djvec_t> jQED_A_bare_0=jQED_qA0_lV0_bare_0+jQED_qAi_lVi_bare_0;
  const vector<djvec_t> jQED_V_bare_0=jQED_qV0_lV0_bare_0+jQED_qVi_lVi_bare_0;

  const vector<djvec_t> jLO_A_bare_ave=(jLO_A_bare_0+jLO_A_bare_1)/2.0;
  const vector<djvec_t> jQED_qA0_lV0_bare_ave=(jQED_qA0_lV0_bare_0+jQED_qA0_lV0_bare_1)/2.0;
  const vector<djvec_t> jQED_qV0_lV0_bare_ave=(jQED_qV0_lV0_bare_0+jQED_qV0_lV0_bare_1)/2.0;
  const vector<djvec_t> jQED_qAi_lVi_bare_ave=(jQED_qAi_lVi_bare_0+jQED_qAi_lVi_bare_1)/2.0;
  const vector<djvec_t> jQED_qVi_lVi_bare_ave=(jQED_qVi_lVi_bare_0+jQED_qVi_lVi_bare_1)/2.0;
  const vector<djvec_t> jQED_A_bare_ave=jQED_qA0_lV0_bare_ave+jQED_qAi_lVi_bare_ave;
  const vector<djvec_t> jQED_V_bare_ave=jQED_qV0_lV0_bare_ave+jQED_qVi_lVi_bare_ave;
  
  size_t ilep=iMLep_of_proc[iproc];
  size_t iQED_mes=iQED_mes_of_proc[iproc];
  
  const double uss=1/(16*sqr(M_PI));
  const double Z3=uss*1.607;
  const double Z4=uss*-3.214;
  const double Wreg_contr_Z3_Z4=Z3-Z4;
  
  //! loop over analysis
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      dbvec_t dA_fr_A_Wreg_dev_1(nens_used);
      dbvec_t dA_fr_A_Wreg_dev_0(nens_used);
      dbvec_t dA_fr_A_Wreg_dev_rdep_1(nens_used);
      dbvec_t dA_fr_A_Wreg_dev_rdep_0(nens_used);
      dbvec_t dA_fr_A_Wreg_dev_Zfact_1(nens_used);
      dbvec_t dA_fr_A_Wreg_dev_Zfact_0(nens_used);
      vector<double> ml_ren(nens_used);
      
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  const size_t iQCD_mes=QED_mes_pars[iQED_mes].iQCD;
	  const size_t ind_QCD=ind_ens_QCD_mes({iens,iQCD_mes});
	  const ens_pars_t &ens=ens_pars[iens];
	  const double aMLep=ens.aMLep[ilep];

	  const size_t tmin=QCD_mes_pars[iQCD_mes].itint;
	  const size_t tmax=ens.T/2-tmin;
	  
	  bi=jack_index[input_an_id][ens.iult];
	  prepare_az(input_an_id);
	  
	  const size_t ib=ens.ib;
	  const size_t ind_proc=ind_ens_proc({iens,iproc});
	  
	  //compute the nasty diagram
	  const dbvec_t QED_qA0_lV0_1=Zv[ib]*dbvec_t(bi,-jQED_qA0_lV0_bare_1[ind_proc]);
	  const dbvec_t QED_qV0_lV0_1=Za[ib]*dbvec_t(bi,jQED_qV0_lV0_bare_1[ind_proc]);
	  const dbvec_t QED_qAi_lVi_1=Zv[ib]*dbvec_t(bi,-jQED_qAi_lVi_bare_1[ind_proc]);
	  const dbvec_t QED_qVi_lVi_1=Za[ib]*dbvec_t(bi,jQED_qVi_lVi_bare_1[ind_proc]);
	  
	  const dbvec_t QED_qA0_lV0_0=Zv[ib]*dbvec_t(bi,-jQED_qA0_lV0_bare_0[ind_proc]);
	  const dbvec_t QED_qV0_lV0_0=Za[ib]*dbvec_t(bi,jQED_qV0_lV0_bare_0[ind_proc]);
	  const dbvec_t QED_qAi_lVi_0=Zv[ib]*dbvec_t(bi,-jQED_qAi_lVi_bare_0[ind_proc]);
	  const dbvec_t QED_qVi_lVi_0=Za[ib]*dbvec_t(bi,jQED_qVi_lVi_bare_0[ind_proc]);
	  
	  const dbvec_t QED_qA0_lV0_ave=Zv[ib]*dbvec_t(bi,-jQED_qA0_lV0_bare_ave[ind_proc]);
	  const dbvec_t QED_qV0_lV0_ave=Za[ib]*dbvec_t(bi,jQED_qV0_lV0_bare_ave[ind_proc]);
	  const dbvec_t QED_qAi_lVi_ave=Zv[ib]*dbvec_t(bi,-jQED_qAi_lVi_bare_ave[ind_proc]);
	  const dbvec_t QED_qVi_lVi_ave=Za[ib]*dbvec_t(bi,jQED_qVi_lVi_bare_ave[ind_proc]);

	  const dbvec_t LO_1=Zv[ib]*dbvec_t(bi,-jLO_A_bare_1[ind_proc]);
	  const dbvec_t LO_0=Zv[ib]*dbvec_t(bi,-jLO_A_bare_0[ind_proc]);
	  const dbvec_t LO_ave=Zv[ib]*dbvec_t(bi,-jLO_A_bare_ave[ind_proc]);
	  const dbvec_t QED_1=Zv[ib]*dbvec_t(bi,-jQED_A_bare_1[ind_proc])+dbvec_t(bi,jQED_V_bare_1[ind_proc])*Za[ib];
	  const dbvec_t QED_0=Zv[ib]*dbvec_t(bi,-jQED_A_bare_0[ind_proc])+dbvec_t(bi,jQED_V_bare_0[ind_proc])*Za[ib];
	  const dbvec_t QED_ave=Zv[ib]*dbvec_t(bi,-jQED_A_bare_ave[ind_proc])+dbvec_t(bi,jQED_V_bare_ave[ind_proc])*Za[ib];
	  
	  const double pi_bare=find_pi(aMLep,ens.aMMes[iQCD_mes]);
	  
	  //cout<<"betal for ensemble "<<ens.path<<": "<<betal.ave_err()<<endl;
	  const dboot_t a=1/lat_par[input_an_id].ainv[ib];
	  const dboot_t L=ens.L*a;
	  const dboot_t MLep=aMLep/a;
	  const dboot_t Mmes=dboot_t(bi,jaM[ind_QCD])/a;
	  
	  //extract dA/A from nasty diagram ratio
	  dbvec_t rat_ext_1=QED_1/LO_1;
	  rat_ext_1[rat_ext_1.size()-1]=rat_ext_1[0]=0.0; //set to zero the contact term
	  dbvec_t rat_ext_0=QED_0/LO_0;
	  rat_ext_0[rat_ext_0.size()-1]=rat_ext_0[0]=0.0;
	  dbvec_t rat_ext_ave=QED_ave/LO_ave;
	  rat_ext_ave[rat_ext_ave.size()-1]=rat_ext_ave[0]=0.0;

	  const dboot_t reno_coeff=Za[ib]/Zv[ib];
	  const dboot_t ZP_fr_ZA=reno_coeff*dboot_t(bi,jZP[ind_QCD])/dboot_t(bi,jZA[ind_QCD]);
	  dboot_t rdep=Wreg_contr_Z3_Z4*ZP_fr_ZA*ratio_Wreg_leptonic_traces(pi_bare,aMLep);
	  
	  const double Zfact=0.75;

	  dboot_t rdep_Zfact=rdep*Zfact;
	  
	  const dbvec_t rat_ext_Wreg_dev_1=rat_ext_1-rat_ext_ave;
	  const dbvec_t rat_ext_Wreg_dev_0=rat_ext_0-rat_ext_ave;
	  
	  const dbvec_t rat_ext_Wreg_dev_rdep_1=rat_ext_Wreg_dev_1-rdep;
	  const dbvec_t rat_ext_Wreg_dev_rdep_0=rat_ext_Wreg_dev_0+rdep;

	  const dbvec_t rat_ext_Wreg_dev_Zfact_1=rat_ext_Wreg_dev_1-rdep_Zfact;
	  const dbvec_t rat_ext_Wreg_dev_Zfact_0=rat_ext_Wreg_dev_0+rdep_Zfact;
	    
	  grace_file_t outA0V0(ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_qA0_lV0_an"+to_string(input_an_id)+".xmg");
	  outA0V0.write_vec_ave_err(QED_qA0_lV0_1.ave_err());
	  outA0V0.write_vec_ave_err(QED_qA0_lV0_0.ave_err());
	  outA0V0.write_vec_ave_err(QED_qA0_lV0_ave.ave_err());
	  grace_file_t outV0V0(ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_qV0_lV0_an"+to_string(input_an_id)+".xmg");
	  outV0V0.write_vec_ave_err(QED_qV0_lV0_1.ave_err());
	  outV0V0.write_vec_ave_err(QED_qV0_lV0_0.ave_err());
	  outV0V0.write_vec_ave_err(QED_qV0_lV0_ave.ave_err());
	  grace_file_t outAiVi(ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_qAi_lVi_an"+to_string(input_an_id)+".xmg");
	  outAiVi.write_vec_ave_err(QED_qAi_lVi_1.ave_err());
	  outAiVi.write_vec_ave_err(QED_qAi_lVi_0.ave_err());
	  outAiVi.write_vec_ave_err(QED_qAi_lVi_ave.ave_err());
	  grace_file_t outViVi(ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_qVi_lVi_an"+to_string(input_an_id)+".xmg");
	  outViVi.write_vec_ave_err(QED_qVi_lVi_1.ave_err());
	  outViVi.write_vec_ave_err(QED_qVi_lVi_0.ave_err());
	  outViVi.write_vec_ave_err(QED_qVi_lVi_ave.ave_err());

	  dA_fr_A_Wreg_dev_1[iens]=constant_fit(rat_ext_Wreg_dev_1,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_1.xmg");
	  dA_fr_A_Wreg_dev_0[iens]=constant_fit(rat_ext_Wreg_dev_0,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_0.xmg");
	  dA_fr_A_Wreg_dev_rdep_1[iens]=constant_fit(rat_ext_Wreg_dev_rdep_1,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_rdep_1.xmg");
	  dA_fr_A_Wreg_dev_rdep_0[iens]=constant_fit(rat_ext_Wreg_dev_rdep_0,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_rdep_0.xmg");
	  dA_fr_A_Wreg_dev_Zfact_1[iens]=constant_fit(rat_ext_Wreg_dev_Zfact_1,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_Zfact_1.xmg");
	  dA_fr_A_Wreg_dev_Zfact_0[iens]=constant_fit(rat_ext_Wreg_dev_Zfact_0,tmin,tmax,ens.path+"/plots_Zfact/proc_"+to_string(iproc)+"_an"+to_string(input_an_id)+"dev_Zfact_0.xmg");
	  
	  ml_ren[iens]=dboot_t(ens.aml*lat_par[input_an_id].ainv[ib]/lat_par[input_an_id].Z[ib]).ave();
	}
      
      //prepare plot
      grace_file_t fit_file("plots_hl/Z_fact_an"+to_string(input_an_id)+".xmg");
      fit_file.set_subtitle("Z_fact");
      fit_file.set_xaxis_label("$$m_{light} (\\overline{MS},2 GeV) [GeV]");
      fit_file.set_yaxis_label("Z_fact");
      fit_file.set_xaxis_max(0.05);
      fit_file.new_data_set();
      
      for(size_t ib=0;ib<nbeta;ib++)
      	{
      	  //make the list of volumes
      	  set<size_t> L_list;
      	  for(size_t idata=0;idata<ens_pars.size();idata++)
      	    if(ens_pars[idata].ib==ib)
      	      L_list.insert(ens_pars[idata].L);
	  
      	  //loop over the list of volumes
      	  for(auto &L : L_list)
      	    {
      	      fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
	      
      	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_1[idata].ave_err());
      	      fit_file.new_data_set();

	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_0[idata].ave_err());
	      fit_file.new_data_set();

	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_rdep_1[idata].ave_err());
	      fit_file.new_data_set();

	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_rdep_0[idata].ave_err());
	      fit_file.new_data_set();

	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_Zfact_1[idata].ave_err());
	      fit_file.new_data_set();

	      for(size_t idata=0;idata<ens_pars.size();idata++)
      	      	if(ens_pars[idata].ib==ib and ens_pars[idata].L==L)
      	      	    fit_file.write_ave_err(dboot_t(ens_pars[idata].aml/lat_par[input_an_id].Z[ib]*lat_par[input_an_id].ainv[ib]).ave(),dA_fr_A_Wreg_dev_Zfact_0[idata].ave_err());
	      fit_file.new_data_set();
      	    }
      	}
    }
}

int main(int narg,char **arg)
{
  int start=time(0);
  
  // cout.precision(16);
  // cout<<zeta(0.27138338825)<<endl;
  // cout<<endl;
  // cout<<zeta(0)<<endl;
  // exit(0);
  
  initialize(narg,arg);
  
  compute_basic_slopes();
  compute_adml_bare();
  
  //preliminary test for factorization
  test_factorization(0);
  
  load_all_hl();
  
  enum{EXCLUDE_IB,INCLUDE_IB};
  
  //loop over process
  dbvec_t tot_corr_proc[nprocess];
  for(size_t iproc=0;iproc<nprocess;iproc++)
    tot_corr_proc[iproc]=compute_corr(iproc,INCLUDE_IB);
  
  const size_t iLEP=0;
  //iFSE_max=0 -> max_ord=1, iFSE_max=1 -> max_ord=2, defined in FSE_max_orders
  // extrapolate_corr(tot_corr_proc[iK]-tot_corr_proc[iPi],{{iK,+1.0},{iPi,-1.0}},STUDY_K_M_PI,iLEP,
  // 		   {hl::FSE::SUB_STDEP::NO,hl::FSE::SUB_STDEP::NO},
  // 		   {hl::FSE::FIT_STDEP::YES,hl::FSE::FIT_STDEP::NO});
  extrapolate_corr(tot_corr_proc[iPi],{{iPi,1.0}},STUDY_PI,iLEP,
  		   {hl::FSE::SUB_STDEP::YES,hl::FSE::SUB_STDEP::NO},
  		   {hl::FSE::FIT_STDEP::YES,hl::FSE::FIT_STDEP::YES});
  // extrapolate_corr(tot_corr_proc[iK],{{iK,1.0}},STUDY_K,iLEP,
  // 		   {hl::FSE::SUB_STDEP::YES,hl::FSE::SUB_STDEP::NO},
  // 		   {hl::FSE::FIT_STDEP::YES,hl::FSE::FIT_STDEP::YES});

  cout<<endl<<"Total time: "<<time(0)-start<<" s"<<endl;
  
  return 0;
}
