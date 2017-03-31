#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <common.hpp>

class qpars_t
{
public:
  size_t iq;
  double dm;
  double eq;
  qpars_t(size_t iq,double dm,double eq) : iq(iq),dm(dm),eq(eq) {}
};
const qpars_t qU(0,1,eu),qD(0,-1,ed),qS(1,0,es),qC(2,0,ec);

//! hold name and consittuent quark for a given meson
class mes_pars_t
{
public:
  size_t iq1;
  size_t iq2;
  double dm1;
  double dm2;
  double eq1;
  double eq2;
  size_t itint;
  string name;
  mes_pars_t(const qpars_t &q1,const qpars_t &q2,size_t itint,const string &name) :
    iq1(q1.iq),iq2(q2.iq),dm1(q1.dm),dm2(q2.dm),eq1(q1.eq),eq2(q2.eq),itint(itint),name(name) {}
};
const vector<mes_pars_t> mes_pars({{qD,qU,0,"PiPlus"},{qS,qU,1,"KPlus"},{qD,qS,1,"K0"},{qD,qC,2,"DPlus"}});
const size_t &nmes=mes_pars.size();
enum{iPiPlus,iKPlus,iK0,iD};

const size_t nmes_tint=3;
class ens_pars_t
{
public:
  size_t iult; //< input in the ultimate file
  size_t ib,T,L;
  double aml;
  int use_for_L;
  string path;
  
  vector<size_t> tmin,tmax;
  ens_pars_t() : tmin(nmes_tint),tmax(nmes_tint) {}
};
vector<ens_pars_t> ens_pars;
size_t nens_used;

const size_t nmass=3,nr=2;
const index_t ind_2pts({{"NMass",nmass},{"NMass",nmass},{"Nr",nr},{"RI",2}});
index_t ind_hl;

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
  
  input.expect({"Ens","beta","aml","L","UseForL","T","TPi","TK","TD","path"});
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

//! compute the critical deltam
djack_t compute_deltam_cr(const ens_pars_t &ens,size_t iq,size_t imes)
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
  
  const size_t itint=mes_pars[imes].itint;
  djack_t deltam_cr=constant_fit(djvec_t(-num_deltam_cr/(2.0*den_deltam_cr)),ens.tmin[itint],ens.tmax[itint],combine("%s/deltam_cr_t.xmg",ens_qpath.c_str()));
  
  return deltam_cr;
}

//! read QED corrections
djvec_t read_QED(const ens_pars_t &ens,size_t imes,const djack_t &deltam_cr,const djvec_t &c_LO)
{
  string ens_qpath=ens.path+"/plots_"+mes_pars[imes].name;
  
  size_t iq1=mes_pars[imes].iq1;
  size_t iq2=mes_pars[imes].iq2;
  double eq1=mes_pars[imes].eq1;
  double eq2=mes_pars[imes].eq2;
  djvec_t c_0T1=read_PP("0T",ens,iq2,iq1,1,RE)*eq1*eq1;
  djvec_t c_0T2=read_PP("0T",ens,iq1,iq2,1,RE)*eq2*eq2;
  djvec_t c_0M1=read_PP("0M",ens,iq2,iq1,1,RE)*eq1*eq1;
  djvec_t c_0M2=read_PP("0M",ens,iq1,iq2,1,RE)*eq2*eq2;
  djvec_t c_LL=read_PP("LL",ens,iq1,iq2,1,RE)*eq1*eq2;
  djvec_t(c_0T1/c_LO).ave_err().write(combine("%s/P5P5_0T1.xmg",ens_qpath.c_str()));
  djvec_t(c_0T2/c_LO).ave_err().write(combine("%s/P5P5_0T2.xmg",ens_qpath.c_str()));
  djvec_t(c_0M1/c_LO).ave_err().write(combine("%s/P5P5_0M1.xmg",ens_qpath.c_str()));
  djvec_t(c_0M2/c_LO).ave_err().write(combine("%s/P5P5_0M2.xmg",ens_qpath.c_str()));
  djvec_t(c_LL/c_LO).ave_err().write(combine("%s/P5P5_LL.xmg",ens_qpath.c_str()));
  djvec_t c=-(c_LL+c_0T1+c_0M1+c_0T2+c_0M2); //minus due to slope definition
  
  djvec_t c_0P1=-(read_PP("0P",ens,iq2,iq1,-1,IM)*eq1*eq1); //remind that i missing means
  djvec_t c_0P2=-(read_PP("0P",ens,iq1,iq2,-1,IM)*eq2*eq2); // to take imag part changed of sign
  djvec_t(c_0P1/c_LO).ave_err().write(combine("%s/P5P5_0P1.xmg",ens_qpath.c_str()));
  djvec_t(c_0P2/c_LO).ave_err().write(combine("%s/P5P5_0P2.xmg",ens_qpath.c_str()));
  djvec_t d=djack_t(-deltam_cr)*(c_0P1+c_0P2);
  
  return djvec_t(c+d)*e2;
}

//! read MASS corrections
djvec_t read_MASS(const ens_pars_t &ens,size_t imes,const djvec_t &c_LO)
{
  string ens_qpath=ens.path+"/plots_"+mes_pars[imes].name;
  
  size_t iq1=mes_pars[imes].iq1;
  size_t iq2=mes_pars[imes].iq2;
  djvec_t c_0S1=read_PP("0S",ens,iq2,iq1,1,RE)*mes_pars[imes].dm1;
  djvec_t c_0S2=read_PP("0S",ens,iq1,iq2,1,RE)*mes_pars[imes].dm2;
  djvec_t(c_0S1/c_LO).ave_err().write(combine("%s/P5P5_0S1.xmg",ens_qpath.c_str()));
  djvec_t(c_0S2/c_LO).ave_err().write(combine("%s/P5P5_0S2.xmg",ens_qpath.c_str()));
  return -(c_0S1+c_0S2);
}

int main(int narg,char **arg)
{
  int start=time(0);
  
  initialize(narg,arg);
  
  cout<<endl<<"Total time: "<<time(0)-start<<" s"<<endl;
  
  index_t ind_ens_mes({{"Ens",nens_used},{"Mes",nmes}});
  size_t nens_mes=ind_ens_mes.max();
  vector<djvec_t> jPP_LO(nens_mes);
  vector<djvec_t> jPP_QED(nens_mes);
  vector<djvec_t> jPP_MASS(nens_mes);
  
  vector<djack_t> ZP(nens_mes);
  vector<djack_t> M(nens_mes);
  vector<djack_t> A_QED(nens_mes);
  vector<djack_t> SL_QED(nens_mes);
  vector<djack_t> A_MASS(nens_mes);
  vector<djack_t> SL_MASS(nens_mes);
  
  //load everything
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_pars_t &ens=ens_pars[iens];
      djack_t deltam_cr=compute_deltam_cr(ens,ilight,iPiPlus);
      size_t TH=ens.T/2;
      
      //load LO for PP
      for(size_t imes=0;imes<nmes;imes++)
	{
	  const size_t ind=ind_ens_mes({iens,imes});
	  const size_t itint=mes_pars[imes].itint;
	  size_t tmin=ens.tmin[itint];
	  size_t tmax=ens.tmax[itint];
	  
	  string plots_path=ens.path+"/plots_"+mes_pars[imes].name;
	  
	  size_t iq1=mes_pars[imes].iq1;
	  size_t iq2=mes_pars[imes].iq2;
	  jPP_LO[ind]=read_PP("00",ens,iq1,iq2,1,RE);
	  jPP_QED[ind]=read_QED(ens,imes,deltam_cr,jPP_LO[ind]);
	  jPP_MASS[ind]=read_MASS(ens,imes,jPP_LO[ind]);
	  
	  two_pts_with_ins_ratio_fit(ZP[ind],M[ind],A_QED[ind],SL_QED[ind],jPP_LO[ind],jPP_QED[ind],TH,tmin,tmax,
				     plots_path+"/effmass_LO1.xmg",plots_path+"/slope_QED.xmg");
	  two_pts_with_ins_ratio_fit(ZP[ind],M[ind],A_MASS[ind],SL_MASS[ind],jPP_LO[ind],jPP_MASS[ind],TH,tmin,tmax,
				     plots_path+"/effmass_LO2.xmg",plots_path+"/slope_MASS.xmg");
	}
    }
  
  index_t ind_ml({{"Input",ninput_an},{"Ens",nens_used}});
  dbvec_t adm_bare(ind_ml.max());
  
  for(size_t ind=0;ind<ind_ml.max();ind++)
    {
      size_t input_an_id=ind_ml(ind)[0];
      size_t iens=ind_ml(ind)[1];
      
      ens_pars_t &ens=ens_pars[iens];
      bi=jack_index[input_an_id][ens.iult];
      size_t ib=ens.ib;
      dboot_t a=1/lat_par[input_an_id].ainv[ib];
      dboot_t ZP=lat_par[input_an_id].Z[ib];
      
      size_t ind_Kplus=ind_ens_mes({iens,iKPlus}),ind_K0=ind_ens_mes({iens,iK0});
      const double phys_dM2K=sqr(MKPLUS)-sqr(MK0);
      dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*ZP*(6.0*log(mu_MS*a)-22.596)/(32.0*sqr(M_PI)));
      
      dboot_t QED_dM2K;
      QED_dM2K=dboot_t(bi,SL_QED[ind_Kplus]-SL_QED[ind_K0]);
      dboot_t ml=ens_pars[iens].aml/ZP/a;
      QED_dM2K+=ml*a*dboot_t(bi,SL_MASS[ind_Kplus]-SL_MASS[ind_K0])/Z_QED;
      QED_dM2K*=2*dboot_t(bi,M[ind_Kplus]);
      QED_dM2K-=dboot_t(bi,FVE_M2(M[ind_Kplus],ens.L));
      QED_dM2K/=sqr(a);
      
      dboot_t QCD_dM2K=phys_dM2K-QED_dM2K;
      
      dboot_t QCD_dM2K_over_adm;
      QCD_dM2K_over_adm=dboot_t(bi,(SL_MASS[ind_Kplus]-SL_MASS[ind_K0])*2*M[ind_Kplus])/sqr(a);
      adm_bare[ind]=QCD_dM2K/QCD_dM2K_over_adm;
      cout<<"dm["<<ind<<"]: "<<dboot_t(adm_bare[ind]/ZP/a).ave_err()<<endl;
  }
  
  return 0;
}
