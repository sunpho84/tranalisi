#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <common.hpp>

const int EVN=1,ODD=-1;//,UNK=0;

typedef int rcombo_2pts_t[2][2];
typedef int rcombo_3pts_t[2][3];

const int SAMER_2PTS=0,OPPOR_2PTS=1;
const rcombo_2pts_t RTAGS_2PTS[2]={{{0,0},{1,1}}, {{0,1},{1,0}}};
const array<string,2> RTAGS_2PTS_NAME={"SAME","OPPO"};

//const int SAMER_SAMER_3PTS=0,SAMER_OPPOR_3PTS=1,OPPOR_SAMER_3PTS=2,OPPOR_OPPOR_3PTS=3;
const array<string,4> RTAGS_3PTS_NAME={"SAME_SAME","SAME_OPPO","OPPO_SAME","OPPO_OPPO"};
const rcombo_3pts_t RTAGS_3PTS[4]={{{0,0,1},{1,1,0}}, {{0,0,0},{1,1,1}}, {{1,0,1},{0,1,0}}, {{0,1,1},{1,0,0}}}; //REV, SPEC, SEQ , remember that SEQ is already reversed

size_t T;
size_t tmin,tmax;

//! return "_var" or ""
string prespaced(const string &var)
{return (var==""?"":"_")+var;}

//! read a single correlator
djvec_t get(const string &name,const int ri,const int tpar)
{
  //cout<<"Reading "<<name<<endl;
  return read_djvec("data/"+name,T,ri).symmetrized(tpar);
}

//! take three vectors and compute derivative
djvec_t derP(const array<djvec_t,3> &corrs,const array<double,3> &coeffs)
{return djvec_t(corrs[1]-corrs[0])/(coeffs[1]-coeffs[0]);}

//! take three vectors and compute derivative
djvec_t derM(const array<djvec_t,3> &corrs,const array<double,3> &coeffs)
{return djvec_t(corrs[2]-corrs[0])/(coeffs[2]-coeffs[0]);}

//! take three vectors and compute derivative
djvec_t der(const array<djvec_t,3> &corrs,const array<double,3> &coeffs)
{
  return
    0.5*djvec_t(corrs[1]-corrs[0])/(coeffs[1]-coeffs[0])+
    0.5*djvec_t(corrs[2]-corrs[0])/(coeffs[2]-coeffs[0]);
}

//! load averaging r
djvec_t load_averaging(const array<string,2> &what,const int ri,const int tpar,const int rpar,const string &diracs,const string &var="")
{
  djvec_t temp[2];
  for(int r=0;r<2;r++)
    temp[r]=get(what[r]+prespaced(var)+"_"+diracs,ri,tpar);
  
  return djvec_t(temp[0]+rpar*temp[1])/(1+abs(rpar));
}

//! load 2pts
djvec_t load_2pts(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const string &var="",const string &suff_rev="",const string &suff_for="")
{
  array<string,2> what;
  const rcombo_2pts_t &rtags=RTAGS_2PTS[irtags];
  for(int r=0;r<2;r++) what[r]="Spect"+to_string(rtags[r][0])+prespaced(suff_rev)+"_Spect"+to_string(rtags[r][1])+prespaced(suff_for);
  
  djvec_t out=load_averaging(what,ri,tpar,rpar,dirtag,var);
  out.ave_err().write("plots/"+dirtag+"_"+RTAGS_2PTS_NAME[irtags]+prespaced(suff_rev)+prespaced(suff_for)+prespaced(var)+".xmg");
  
  return out;
}

//! load 2pts averaging on space variable 1 2 3
djvec_t load_2pts_spave(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const string &var="",const string &suff_rev="",const string &suff_for="")
{
  djvec_t out(T/2+1);
  
  out=0.0;
  for(char c : {'1','2','3'})
    out+=load_2pts(combine(dirtag.c_str(),c,c),irtags,ri,tpar,rpar,var,combine(suff_rev.c_str(),c),combine(suff_for.c_str(),c));
  out/=3.0;
  
  out.ave_err().write("plots/"+combine(dirtag.c_str(),'K','K')+"_"+RTAGS_2PTS_NAME[irtags]+
		      prespaced(combine(suff_rev.c_str(),'K'))+
		      prespaced(combine(suff_for.c_str(),'K'))+
		      prespaced(var)+".xmg");
  
  return out;
}

//! compute the derivative od 2pts correlators w.r.t a given variation
djvec_t der_2pts(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const array<double,3> &coeffs,const string &var,const string &suff_rev="",const string &suff_for="",djvec_t(*load_fun)(const string&,const int,const int,const int,const int,const string&,const string&,const string&)=load_2pts)
{
  array<string,3> vars={"",var+"P",var+"M"};
  array<djvec_t,3> corrs;
  for(int i=0;i<3;i++) corrs[i]=load_fun(dirtag,irtags,ri,tpar,rpar,vars[i],suff_rev,suff_for);
  
  derP(corrs,coeffs).ave_err().write("plots/"+combine(dirtag.c_str(),'K','K')+"_"+RTAGS_2PTS_NAME[irtags]+
				     prespaced(combine(suff_rev.c_str(),'K'))+
				     prespaced(combine(suff_for.c_str(),'K'))+
				     "_der_"+var+"P.xmg");
  
  derM(corrs,coeffs).ave_err().write("plots/"+combine(dirtag.c_str(),'K','K')+"_"+RTAGS_2PTS_NAME[irtags]+
				     prespaced(combine(suff_rev.c_str(),'K'))+
				     prespaced(combine(suff_for.c_str(),'K'))+
				     "_der_"+var+"M.xmg");
  
  djvec_t out=der(corrs,coeffs);
  out.ave_err().write("plots/"+combine(dirtag.c_str(),'K','K')+"_"+RTAGS_2PTS_NAME[irtags]+
		      prespaced(combine(suff_rev.c_str(),'K'))+
		      prespaced(combine(suff_for.c_str(),'K'))+
		      "_der_"+var+".xmg");
  
  return out;
}

djvec_t der_2pts_spave(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const array<double,3> &coeffs,const string &var,const string &suff_rev="",const string &suff_for="")
{return der_2pts(dirtag,irtags,ri,tpar,rpar,coeffs,var,suff_rev,suff_for,load_2pts_spave);}

//! load 3pts (V0P5)
djvec_t load_3pts(const int irtags,const int tpar,const int rpar,const string &var="")
{
  array<string,2> what;
  const rcombo_3pts_t &rtags=RTAGS_3PTS[irtags];
  
  for(int r=0;r<2;r++) what[r]="Spect"+to_string(rtags[r][0])+"_Seq"+to_string(rtags[r][2])+to_string(rtags[r][1]); //note the r order
  
  djvec_t out=load_averaging(what,RE,tpar,rpar,"V0P5",var);
  out.ave_err().write("plots/P5V0P5_"+RTAGS_3PTS_NAME[irtags]+prespaced(var)+".xmg");
  
  return out;
}

//! compute the derivative of 3pts correlators w.r.t a given variation
djvec_t der_3pts(const int irtags,const int tpar,const int rpar,const array<double,3> &coeffs,const string &var)
{
  array<string,3> vars={"",var+"P",var+"M"};
  array<djvec_t,3> corrs;
  for(int i=0;i<3;i++) corrs[i]=load_3pts(irtags,tpar,rpar,vars[i]);
  
  djvec_t out=der(corrs,coeffs);
  out.ave_err().write("plots/P5V0P5_"+RTAGS_3PTS_NAME[irtags]+"_der_"+var+".xmg");
  
  return out;
}

//! return the three pts behaviour
djack_t three_pts_dt_reco(const djack_t &Zso,const djack_t &Mso,const djack_t &Zsi,const djack_t &Msi,int t)
{return Zso*Zsi/(4*Mso*Msi)*exp(-Mso*t)*exp(-Msi*(T/2.0-t))*(Mso+Msi);}

//! compute the Z from 2pts decay constant
void z_from_2pts_dec(djack_t &za,djack_t &zv,const djvec_t &P5P5_SAME,const djvec_t &A0P5_SAME,const djvec_t &P5P5_OPPO,const djvec_t &A0P5_OPPO,const double m,const string &tag)
{
  djack_t ZP5P5_SAME,MP5P5_SAME;
  djack_t ZP5P5_OPPO,MP5P5_OPPO;
  djack_t ZA0P5_SAME,MA0P5_SAME;
  djack_t ZA0P5_OPPO,MA0P5_OPPO;
  
  two_pts_fit(ZP5P5_SAME,MP5P5_SAME,P5P5_SAME,T/2,tmin,tmax,"plots/effmass_"+tag+"_P5P5_same.xmg");
  two_pts_fit(ZP5P5_OPPO,MP5P5_OPPO,P5P5_OPPO,T/2,tmin,tmax,"plots/effmass_"+tag+"_P5P5_oppo.xmg");
  two_pts_fit(ZA0P5_SAME,MA0P5_SAME,A0P5_SAME,T/2,tmin,tmax,"plots/effmass_"+tag+"_A0P5_same.xmg","",-1);
  two_pts_fit(ZA0P5_OPPO,MA0P5_OPPO,A0P5_OPPO,T/2,tmin,tmax,"plots/effmass_"+tag+"_A0P5_oppo.xmg","",-1);
  
  djack_t ZP_SAME=sqrt(ZP5P5_SAME);
  djack_t ZA_SAME=ZA0P5_SAME/ZP_SAME;
  djack_t f_SAME=ZP_SAME*(2*m)/sqr(MP5P5_SAME);
  djack_t f_SAME_bare=ZA_SAME/MP5P5_SAME;
  
  djack_t ZA_OPPO=ZA0P5_OPPO/sqrt(ZP5P5_OPPO);
  djack_t f_OPPO_bare=ZA_OPPO/MP5P5_OPPO;
  
  za=f_SAME/f_OPPO_bare;
  zv=f_SAME/f_SAME_bare;
}

//! perform a fit to determine the slope and mass of 2 2pts and 2 ins
template <class TV,class TS=typename TV::base_type> void two_two_pts_with_ins_ratio_fit(TV &Z2,TS &M,TV &DZ2_fr_Z2,TS &SL,const vector<TV> &corr,const vector<TV> &corr_ins,size_t TH,size_t tmin,size_t tmax,vector<string> path={},vector<string> path_ins={},vector<int> par={})
{
  size_t ncorrs=corr.size();
  if(corr.size()!=corr_ins.size()) CRASH("Sizes of corr %zu and corr_ins %zu not matching",corr.size(),corr_ins.size());
  
  par.resize(ncorrs,1);
  path.resize(ncorrs,"");
  path_ins.resize(ncorrs,"");
  
  //perform a preliminary fit
  bool real_fit=true;
  TV eff_mass=effective_mass(corr[0],TH,par[0]);
  M=constant_fit(eff_mass,tmin,tmax,"/tmp/test_mass.xmg");
  TV eff_slope=effective_slope(TV(corr_ins[0]/corr[0]),eff_mass,TH,par[0]);
  SL=constant_fit(eff_slope,tmin,tmax,"/tmp/test_slope.xmg");
  real_fit&=(SL.err()!=0.0);
  
  for(size_t icorr=0;icorr<ncorrs;icorr++)
    {
      TV eff_sq_coupling=effective_squared_coupling(corr[icorr],eff_mass,TH,par[icorr]);
      Z2[icorr]=constant_fit(eff_sq_coupling,tmin,tmax,"/tmp/test_sq_coupling"+to_string(icorr)+".xmg");
      TV eff_slope_offset=effective_squared_coupling_rel_corr(TV(corr_ins[icorr]/corr[icorr]),eff_mass,eff_slope,TH,par[icorr]);
      DZ2_fr_Z2[icorr]=constant_fit(eff_slope_offset,tmin,tmax,"/tmp/test_sq_coupling_rel_corr"+to_string(icorr)+".xmg");
      real_fit&=(DZ2_fr_Z2[icorr].ave()!=0.0);
    }
  
  if(real_fit)
    {
      //! fit for real
      size_t iel=0;
      
      minimizer_pars_t pars;
      size_t iM=0;pars.add("M",M.ave(),M.err());
      size_t iSl=1;pars.add("SL",SL.ave(),SL.err());
      vector<size_t> iZ2(ncorrs),iDZ2_fr_Z2(ncorrs);
      for(size_t icorr=0;icorr<ncorrs;icorr++) {iZ2[icorr]=2+icorr;pars.add("Z2_"+to_string(icorr),Z2[icorr].ave(),Z2[icorr].err());}
      for(size_t icorr=0;icorr<ncorrs;icorr++) {iDZ2_fr_Z2[icorr]=2+icorr+ncorrs;pars.add("DZ2_fr_Z2_"+to_string(icorr),DZ2_fr_Z2[icorr].ave(),DZ2_fr_Z2[icorr].err());}
      
      auto x=vector_up_to<double>(corr[0].size());
      vector<function<double(const vector<double> &p,double x)>> funs;
      vector<djvec_t> corrs;
      for(size_t icorr=0;icorr<ncorrs;icorr++)
	{
	  funs.push_back(two_pts_corr_fun_t(TH,par[icorr]));
	  corrs.push_back(corr[icorr]);
	}
      for(size_t icorr=0;icorr<ncorrs;icorr++)
	{
	  funs.push_back(two_pts_corr_with_ins_fun_t(TH,par[icorr]));
	  corrs.push_back(corr_ins[icorr]/corr[icorr]);
	}
      multi_ch2_t<TV> two_pts_fit_obj(vector<vector<double>>(2*ncorrs,x),
				      vector<size_t>(2*ncorrs,tmin),
				      vector<size_t>(2*ncorrs,tmax),
				      corrs,
				      funs,
				      [ncorrs,iM,&iZ2,&iDZ2_fr_Z2,iSl]
				      (const vector<double> &p,size_t icorr)
				      {
					if(icorr<ncorrs) return vector<double>({p[iZ2[icorr%ncorrs]],p[iM]});
					else             return vector<double>({p[iM],p[iDZ2_fr_Z2[icorr%ncorrs]],p[iSl]});
				      }
				      ,iel);
      
      //parameters to fit
      minimizer_t minimizer(two_pts_fit_obj,pars);
      
      for(iel=0;iel<corr[0][0].size();iel++)
	{
	  //minimize and print the result
	  vector<double> par_min=minimizer.minimize();
	  M[iel]=par_min[iM];
	  SL[iel]=par_min[iSl];
	  for(size_t icorr=0;icorr<ncorrs;icorr++)
	    {
	      Z2[icorr][iel]=par_min[iZ2[icorr]];
	      DZ2_fr_Z2[icorr][iel]=par_min[iDZ2_fr_Z2[icorr]];
	    }
	}
      
      //write plots
      for(size_t icorr=0;icorr<ncorrs;icorr++)
	{
	  if(path[icorr]!="") write_constant_fit_plot(path[icorr],tmin,tmax,M,effective_mass(corr[icorr],TH,par[icorr]));
	  if(path_ins[icorr]!="") write_fit_plot(path_ins[icorr],tmin,tmax,[M,&DZ2_fr_Z2,SL,TH,icorr,&par](double x)->TS
				      {return two_pts_corr_with_ins_ratio_fun(M,DZ2_fr_Z2[icorr],SL,TH,x,par[icorr]);},TV(corr_ins[icorr]/corr[icorr]));
	}
    }
  else CRASH("cannot fit");
}

int main(int narg,char **arg)
{
  string name="input.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  T=input.read<size_t>("T");
  tmin=input.read<size_t>("tmin_2pts");
  tmax=input.read<size_t>("tmax_2pts");
  
  size_t temp_njacks=input.read<size_t>("NJacks");
  set_njacks(temp_njacks);
  
  bool ant=input.read<bool>("ant");
  double use_mass=input.read<double>("use_mass");
  
  const double k=input.read<double>("k");
  double kp=0,km=0;
  if(ant)
    {
      kp=input.read<double>("kp");
      km=input.read<double>("km");
    }
  
  const double ma=input.read<double>("m");
  
  double map=0,mam=0;
  if(ant)
    {
      map=input.read<double>("mp");
      mam=input.read<double>("mm");
    }
  const size_t ib=input.read<size_t>("ib");
  
  const size_t WI_2pts=input.read<size_t>("WI_2pts");
  const size_t WI_3pts=0;

  const double charge=input.read<double>("charge");

  const double qf2_der=sqr(charge/3),qf2_phys=1.0/9.0;
  const array<double,3> e2_var={0,qf2_der*e2,qf2_der*e2};
  const double e2_phys=e2*qf2_phys;
  
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  const array<double,3> ka_var={mc,mcp,mcm};
  array<double,3> ma_var={ma,map,mam};
  
  const double a[3]={0.45063,0.414834,0.31476};
  const double pred_Zp_coeff=6.0*log(mu_MS*a[ib])-22.5954;
  const double pred_Zs_coeff=6.0*log(mu_MS*a[ib])-12.9524;
  const double pred_Zv_coeff=-20.6178;
  const double pred_Za_coeff=-15.7963;
  const double pred_Zp_m_Zs_coeff=pred_Zp_coeff-pred_Zs_coeff;
  const double pred_Za_m_Zv_coeff=pred_Za_coeff-pred_Zv_coeff;
  //const double pred_Zp=pred_Zp_coeff*e2/sqr(4*M_PI);
  //const double pred_Zs=pred_Zs_coeff*e2/sqr(4*M_PI);
  const double pred_Zv=pred_Zv_coeff*e2/sqr(4*M_PI);
  const double pred_Za=pred_Za_coeff*e2/sqr(4*M_PI);
  
  djvec_t V0P5_SAME=load_2pts("V0P5",SAMER_2PTS,IM,ODD,ODD);
  djvec_t V0P5_OPPO=load_2pts("V0P5",OPPOR_2PTS,IM,ODD,ODD);
  V0P5_SAME.ave_err().write("plots/V0P5_SAME.xmg");
  V0P5_OPPO.ave_err().write("plots/V0P5_OPPO.xmg");
  djvec_t P5P5_SAME=load_2pts("P5P5",SAMER_2PTS,RE,EVN,EVN);
  djvec_t P5P5_OPPO=load_2pts("P5P5",OPPOR_2PTS,RE,EVN,EVN);
  djvec_t mcr_SAME=backward_derivative(V0P5_SAME)/P5P5_SAME;
  djvec_t mcr_OPPO=backward_derivative(V0P5_OPPO)/P5P5_OPPO;
  mcr_SAME.ave_err().write("plots/mcr_SAME.xmg");
  mcr_OPPO.ave_err().write("plots/mcr_OPPO.xmg");
  
  //compute deltamcr
  djvec_t E2_V0P5_SAME;
  if(ant) E2_V0P5_SAME=der_2pts("V0P5",SAMER_2PTS,IM,ODD,ODD,e2_var,"E2E2");
  else    E2_V0P5_SAME=
	    load_2pts("V0P5_M0",SAMER_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_0M",SAMER_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_T0",SAMER_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_0T",SAMER_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_LL",SAMER_2PTS,IM,ODD,ODD);
  E2_V0P5_SAME=backward_derivative(E2_V0P5_SAME);
  E2_V0P5_SAME.ave_err().write("plots/V0P5_SAME_E2.xmg");
  
  djvec_t KA_V0P5_SAME;
  if(ant) KA_V0P5_SAME=der_2pts("V0P5",SAMER_2PTS,IM,ODD,ODD,ka_var,"KAKA");
  else    KA_V0P5_SAME=load_2pts("V0P5_P0",SAMER_2PTS,RE,ODD,EVN)-load_2pts("V0P5_0P",SAMER_2PTS,RE,ODD,EVN);
  KA_V0P5_SAME=backward_derivative(KA_V0P5_SAME);
  KA_V0P5_SAME.ave_err().write("plots/V0P5_SAME_KA.xmg");
  
  //compute deltamcr
  djvec_t E2_V0P5_OPPO;
  if(ant) E2_V0P5_OPPO=der_2pts("V0P5",OPPOR_2PTS,IM,ODD,ODD,e2_var,"E2E2");
  else    E2_V0P5_OPPO=
	    load_2pts("V0P5_M0",OPPOR_2PTS,IM,ODD,ODD)-
	    load_2pts("V0P5_0M",OPPOR_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_T0",OPPOR_2PTS,IM,ODD,ODD)-
	    load_2pts("V0P5_0T",OPPOR_2PTS,IM,ODD,ODD)+
	    load_2pts("V0P5_LL",OPPOR_2PTS,IM,ODD,ODD);
  E2_V0P5_OPPO.ave_err().write("plots/V0P5_OPPO_E2.xmg");
  
  djvec_t KA_V0P5_OPPO;
  if(ant) KA_V0P5_OPPO=der_2pts("V0P5",OPPOR_2PTS,IM,ODD,ODD,ka_var,"KAKA");
  else    KA_V0P5_OPPO=load_2pts("V0P5_P0",OPPOR_2PTS,RE,ODD,EVN)-load_2pts("V0P5_0P",OPPOR_2PTS,RE,ODD,EVN);
  KA_V0P5_OPPO.ave_err().write("plots/V0P5_OPPO_KA.xmg");
  
  djvec_t deltam_cr_same_t=E2_V0P5_SAME/KA_V0P5_SAME;
  size_t tmax_dmcr[3]={20,28,40};
  djack_t deltam_cr_same=constant_fit(deltam_cr_same_t,tmin,tmax_dmcr[ib],"plots/deltam_cr_same.xmg");
  cout<<"deltam_cr_same: "<<deltam_cr_same<<endl;
  
  djvec_t deltam_cr_oppo_t=E2_V0P5_OPPO/KA_V0P5_OPPO;
  djack_t deltam_cr_oppo=constant_fit(deltam_cr_oppo_t,tmin,tmax_dmcr[ib],"plots/deltam_cr_oppo.xmg");
  cout<<"deltam_cr_oppo: "<<deltam_cr_oppo<<endl;

  djack_t deltam_cr_gm2;
  string deltam_cr_gm2_path="deltam_cr_light.bin";
  if(file_exists(deltam_cr_gm2_path))
    {
      deltam_cr_gm2.bin_read(deltam_cr_gm2_path);
      cout<<"deltam_cr_gm2: "<<deltam_cr_gm2<<endl;
    }
  
  double deltam=ma/sqr(4*M_PI)*pred_Zp_coeff;
  djack_t deltam_cr=deltam_cr_oppo;
  
  if(WI_3pts)
    {
      djvec_t Z(2),Z2(2),M(2);
      
      djvec_t Z_CORR(2),Z2_CORR(2),M_CORR(2);
      djvec_t LO_P5P5[2],CORR_P5P5[2];
      for(int rdiff=0;rdiff<2;rdiff++)
	{
	  LO_P5P5[rdiff]=load_2pts("P5P5",rdiff,RE,EVN,EVN);
	  two_pts_fit(Z2[rdiff],M[rdiff],LO_P5P5[rdiff],T/2,tmin,tmax,"plots/effmass_P5P5_"+RTAGS_2PTS_NAME[rdiff]+".xmg");
	  Z[rdiff]=sqrt(Z2[rdiff]);
	  
	  djvec_t E2_P5P5=der_2pts("P5P5",rdiff,RE,EVN,EVN,e2_var,"E2");
	  djvec_t KA_P5P5=der_2pts("P5P5",rdiff,RE,EVN,EVN,ka_var,"KA");
	  
	  CORR_P5P5[rdiff]=LO_P5P5[rdiff]+e2_phys*djvec_t(E2_P5P5-deltam_cr*KA_P5P5);
	  CORR_P5P5[rdiff].ave_err().write("plots/P5P5_"+RTAGS_2PTS_NAME[rdiff]+"_corrected.xmg");
	  two_pts_fit(Z2_CORR[rdiff],M_CORR[rdiff],CORR_P5P5[rdiff],T/2,tmin,tmax,"plots/effmass_P5P5_"+RTAGS_2PTS_NAME[rdiff]+"_corrected.xmg");
	  Z_CORR[rdiff]=sqrt(Z2_CORR[rdiff]);
	  
	  cout<<RTAGS_2PTS_NAME[rdiff]<<endl;
	  cout<<" Mass: "<<M[rdiff].ave_err()<<endl;
	  cout<<" Mass corr: "<<M_CORR[rdiff].ave_err()<<endl;
	  cout<<" Mass diff: "<<djack_t(M_CORR[rdiff]-M[rdiff]).ave_err()<<endl;
	}
      
      for(int rdiff_so=0;rdiff_so<2;rdiff_so++)
	for(int rdiff_si=0;rdiff_si<2;rdiff_si++)
	  {
	    int rdiff_tot=rdiff_si+2*rdiff_so;
	    cout<<"Case: "<<RTAGS_3PTS_NAME[rdiff_tot]<<" "<<((rdiff_so==rdiff_si)?"Zv":"Za")<<endl;
	    
	    djvec_t LO_3=load_3pts(rdiff_tot,ODD,EVN);
	    
	    //LO
	    djvec_t Zren_t;
	    djack_t Zren;
	    {
	      //reconstruct 3pts dt
	      djvec_t LO_3_dt_an(T/2+1),LO_3_dt_nu(T/2+1);
	      for(size_t t=0;t<=T/2;t++)
		{
		  LO_3_dt_an[t]=three_pts_dt_reco(Z[rdiff_so],M[rdiff_so],Z[rdiff_si],M[rdiff_si],t);
		  LO_3_dt_nu[t]=(M[rdiff_so]+M[rdiff_si])*
		    (LO_P5P5[rdiff_so][t]-two_pts_corr_fun(Z2[rdiff_so],M[rdiff_so],T/2.0,T-t,0))*
		    (LO_P5P5[rdiff_si][T/2-t]-two_pts_corr_fun(Z2[rdiff_si],M[rdiff_si],T/2.0,T/2.0+t,0))
		    /(Z[rdiff_so]*Z[rdiff_si]);
		}
	      //cout<<"Check "<<RTAGS_3PTS_NAME[rdiff_tot]<<"reco, "<<LO_3_dt[0].ave_err()<<" "<<djack_t(LO_P5P5[rdiff_so][T/2]/2).ave_err()<<endl;
	      LO_3_dt_an.ave_err().write("plots/LO_3_dt_an_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	      LO_3_dt_nu.ave_err().write("plots/LO_3_dt_nu_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	      
	      // cout<<"Analytic: "<<corr(djack_t(-LO_3_dt[T/4]),LO_3[T/4])<<endl;
	      // cout<<"Numeric: "<<corr(djack_t(-LO_P5P5[rdiff_so][T/2]/2),LO_3[T/4])<<endl;
	      
	      Zren_t=-LO_3_dt_an/LO_3;
	      Zren_t.ave_err().write("plots/Z_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	      Zren=Zren_t[T/4];
	      
	      cout<<" ZWI: "<<smart_print(Zren.ave_err())<<endl;
	      cout<<" ZRI: "<<smart_print(((rdiff_so==rdiff_si)?Zv_ae:Za_ae)[0][ib])<<endl;
	    }
	    
	    //QED
	    djvec_t CORR_Zren_t;
	    djack_t CORR_Zren;
	    {
	      //build corrected 3pts
	      djvec_t E2_3=der_3pts(rdiff_tot,ODD,EVN,e2_var,"E2");
	      djvec_t KA_3=der_3pts(rdiff_tot,ODD,EVN,ka_var,"KA");
	      djvec_t CORR_3=LO_3+e2_phys*djvec_t(E2_3-deltam_cr*KA_3);
	      
	      //write the aperiodic effective mass of 3pts
	      {
		grace_file_t eff_file("plots/effmass_3pts_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
		eff_file.write_vec_ave_err(effective_mass(LO_3,T/2,0,1).ave_err());
		eff_file.write_constant_band(0,T/2.0,djack_t(M[rdiff_so]-M[rdiff_si]));
	      }
	      
	      //write the aperiodic effective mass of 3pts corrected minus non corrected
	      {
		grace_file_t eff_corr_file("plots/effmass_3pts_corr_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
		eff_corr_file.write_vec_ave_err(djvec_t(effective_mass(CORR_3,T/2,0,1)-effective_mass(LO_3,T/2,0,1)).ave_err());
		eff_corr_file.write_constant_band(0,T/2.0,djack_t(M_CORR[rdiff_so]-M_CORR[rdiff_si]-M[rdiff_so]+M[rdiff_si]));
	      }
	      
	      //reconstruct 3pts dt
	      djvec_t CORR_3_dt(T/2+1);
	      for(size_t t=0;t<=T/2;t++) CORR_3_dt[t]=three_pts_dt_reco(Z_CORR[rdiff_so],M_CORR[rdiff_so],Z_CORR[rdiff_si],M_CORR[rdiff_si],t);
	      //cout<<"Check CORR "<<RTAGS_3PTS_NAME[rdiff_tot]<<"reco, "<<CORR_3_dt[0].ave_err()<<" "<<djack_t(CORR_P5P5[rdiff_so][T/2]/2).ave_err()<<endl;
	      CORR_3_dt.ave_err().write("plots/CORR_3_dt_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	      
	      CORR_Zren_t=-CORR_3_dt/CORR_3;
	      CORR_Zren=constant_fit(CORR_Zren_t,T/4-2,T/4+2,"plots/CORR_Z_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	    }
	    
	    double ZQED_pred=(rdiff_so==rdiff_si)?pred_Zv:pred_Za;
	    djvec_t ZQED_fact_t=(CORR_Zren_t-Zren_t)/(Zren_t*qf2_phys);
	    djack_t ZQED_fact=constant_fit(ZQED_fact_t,T/4-2,T/4+2,"plots/ZQED_fact_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	    
	    cout<<" ZQED assuming factoriz: "<<smart_print(ZQED_fact.ave_err())<<endl;
	    cout<<" ZQED according to AOKI: "<<ZQED_pred<<endl;
	    
	    cout<<" Correction to factorization hypotesis: "<<smart_print(djack_t(ZQED_fact/ZQED_pred-1.0).ave_err())<<endl;
	  }
    }
  
  if(WI_2pts)
    {
      int RSAME=0,ROPPO=1;
      
      //determination of Za from 2pts WI
      {
	cout<<"Determining Za and Zav from 2pts decay constant"<<endl;
	
	auto load_LO_CORR=[&e2_var,&ka_var,&ma_var,&e2_phys,&deltam_cr,&deltam,&use_mass,&ant]
	  (djvec_t &LO,djvec_t &DELTA,djvec_t &CORR,const string &name,const int &rcombo,const int &reim,const int &tpar,const int &rpar)
	  {
	    LO=load_2pts(name,rcombo,reim,tpar,rpar);
	    LO.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_LO.xmg");
	    djvec_t E2;
	    if(ant) E2=der_2pts(name,rcombo,reim,tpar,rpar,e2_var,"E2E2");
	    else    E2=
		      load_2pts(name+"_LL",rcombo,reim,tpar,rpar)+
		      load_2pts(name+"_0M",rcombo,reim,tpar,rpar)+
		      load_2pts(name+"_M0",rcombo,reim,tpar,rpar)+
		      load_2pts(name+"_0T",rcombo,reim,tpar,rpar)+
		      load_2pts(name+"_T0",rcombo,reim,tpar,rpar);
	    E2.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_E2.xmg");
	    
	    djvec_t KA;
	    if(ant) KA=der_2pts(name,rcombo,reim,tpar,rpar,ka_var,"KAKA");
	    else
	      if(rcombo==0) KA=load_2pts(name+"_0P",rcombo,!reim,tpar,!rpar)-load_2pts(name+"_P0",rcombo,!reim,tpar,!rpar);
	      else          KA=-load_2pts(name+"_0P",rcombo,!reim,tpar,!rpar)-load_2pts(name+"_P0",rcombo,!reim,tpar,!rpar);
	    KA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_KA.xmg");
	    
	    djvec_t MA;
	    if(ant) MA=der_2pts(name,rcombo,reim,tpar,rpar,ma_var,"MAMA");
	    else    MA=-(load_2pts(name+"_0S",rcombo,reim,tpar,rpar)+load_2pts(name+"_S0",rcombo,reim,tpar,rpar));
	    MA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_MA.xmg");
	    
	    DELTA=E2-deltam_cr*KA+deltam*MA*use_mass;
	    DELTA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_DELTA.xmg");
	    CORR=LO+e2_phys*DELTA;
	  };
	
	djvec_t LO_P5P5_SAME,DELTA_P5P5_SAME,CORR_P5P5_SAME;
	load_LO_CORR(LO_P5P5_SAME,DELTA_P5P5_SAME,CORR_P5P5_SAME,"P5P5",RSAME,RE,EVN,EVN);
	djvec_t LO_P5P5_OPPO,DELTA_P5P5_OPPO,CORR_P5P5_OPPO;
	load_LO_CORR(LO_P5P5_OPPO,DELTA_P5P5_OPPO,CORR_P5P5_OPPO,"P5P5",ROPPO,RE,EVN,EVN);
	djvec_t LO_A0P5_SAME,DELTA_A0P5_SAME,CORR_A0P5_SAME;
	load_LO_CORR(LO_A0P5_SAME,DELTA_A0P5_SAME,CORR_A0P5_SAME,"A0P5",RSAME,RE,ODD,EVN);
	djvec_t LO_A0P5_OPPO,DELTA_A0P5_OPPO,CORR_A0P5_OPPO;
	load_LO_CORR(LO_A0P5_OPPO,DELTA_A0P5_OPPO,CORR_A0P5_OPPO,"A0P5",ROPPO,RE,ODD,EVN);
	
       //global fit
	
	djack_t M_OS,SL_OS;
	djvec_t Z2_OS(2),DZ2_fr_Z2_OS(2);
	two_two_pts_with_ins_ratio_fit(Z2_OS,M_OS,DZ2_fr_Z2_OS,SL_OS,vector<djvec_t>({LO_P5P5_OPPO,LO_A0P5_OPPO}),vector<djvec_t>({DELTA_P5P5_OPPO,DELTA_A0P5_OPPO}),T/2,tmin,tmax,{"plots/P5P5_OPPO_effmass.xmg","plots/A0P5_OPPO_effmass.xmg"},{"plots/P5P5_OPPO_slope.xmg","plots/A0P5_OPPO_slope.xmg"},{1,-1});
	djack_t ZP_OS=sqrt(Z2_OS[0]);
	djack_t ZA_OS=Z2_OS[1]/ZP_OS;
	djack_t DZP_fr_ZP_OS=DZ2_fr_Z2_OS[0]/2;
	djack_t DZA_fr_ZA_OS=DZ2_fr_Z2_OS[1]-DZP_fr_ZP_OS;
	
	djack_t M_TM,SL_TM;
	djvec_t Z2_TM(2),DZ2_fr_Z2_TM(2);
	two_two_pts_with_ins_ratio_fit(Z2_TM,M_TM,DZ2_fr_Z2_TM,SL_TM,vector<djvec_t>({LO_P5P5_SAME,LO_A0P5_SAME}),vector<djvec_t>({DELTA_P5P5_SAME,DELTA_A0P5_SAME}),T/2,tmin,tmax,{"plots/P5P5_SAME_effmass.xmg","plots/A0P5_SAME_effmass.xmg"},{"plots/P5P5_SAME_slope.xmg","plots/A0P5_SAME_slope.xmg"},{1,-1});
	djack_t ZP_TM=sqrt(Z2_TM[0]);
	djack_t ZA_TM=Z2_TM[1]/ZP_TM;
	djack_t DZP_fr_ZP_TM=DZ2_fr_Z2_TM[0]/2;
	djack_t DZA_fr_ZA_TM=DZ2_fr_Z2_TM[1]-DZP_fr_ZP_TM;
	
	// cout<<"Check consistency TM M,  P: "<<M_PP_TM<<", A: "<<M_AP_TM<<", mixed: "<<M_TM<<endl;
	// cout<<"Check consistency OS M,  P: "<<M_PP_OS<<", A: "<<M_AP_OS<<", mixed: "<<M_OS<<endl;
	
	//djack_t za_QED_fact_exp_from_dec=DZP_fr_ZP_TM+2*SL_TM/M_TM-SL_OS/M_OS-DZA_fr_ZA_OS;
	cout<<"SL_OS: "<<SL_OS<<", SL_TM: "<<SL_TM<<endl;
	cout<<"M_OS: "<<M_OS<<", M_TM: "<<M_TM<<endl;
	cout<<"SL/M_OS: "<<djack_t(SL_OS/M_OS)<<", SL/M_TM: "<<djack_t(SL_TM/M_TM)<<endl;
	cout<<"DP_OS: "<<DZP_fr_ZP_OS<<", DP_TM: "<<DZP_fr_ZP_TM<<endl;
	cout<<"DA0_OS: "<<DZA_fr_ZA_OS<<", DA0_TM: "<<DZA_fr_ZA_TM<<endl;
	
	if(0)
	{
	two_pts_with_ins_ratio_fit(Z2_OS[0],M_OS,DZ2_fr_Z2_OS[0],SL_OS,LO_P5P5_OPPO,DELTA_P5P5_OPPO,T/2,tmin,tmax,"plots/qP5P5_OPPO_effmass.xmg","plots/qP5P5_OPPO_slope.xmg");
	two_pts_with_ins_ratio_fit(Z2_TM[0],M_TM,DZ2_fr_Z2_TM[0],SL_TM,LO_P5P5_SAME,DELTA_P5P5_SAME,T/2,tmin,tmax,"plots/qP5P5_SAME_effmass.xmg","plots/qP5P5_SAME_slope.xmg");
	two_pts_with_ins_ratio_fit(Z2_OS[1],M_OS,DZ2_fr_Z2_OS[1],SL_OS,LO_A0P5_OPPO,DELTA_A0P5_OPPO,T/2,tmin,tmax,"plots/qA0P5_OPPO_effmass.xmg","plots/qA0P5_OPPO_slope.xmg",-1);
	two_pts_with_ins_ratio_fit(Z2_TM[1],M_TM,DZ2_fr_Z2_TM[1],SL_TM,LO_A0P5_SAME,DELTA_A0P5_SAME,T/2,tmin,tmax,"plots/qA0P5_SAME_effmass.xmg","plots/qA0P5_SAME_slope.xmg",-1);

	ZP_TM=sqrt(Z2_TM[0]);
	ZA_TM=Z2_TM[1]/ZP_TM;
	DZP_fr_ZP_TM=DZ2_fr_Z2_TM[0]/2;
	DZA_fr_ZA_TM=DZ2_fr_Z2_TM[1]-DZP_fr_ZP_TM;
	
	cout<<"SL_OS: "<<SL_OS<<", SL_TM: "<<SL_TM<<endl;
	cout<<"M_OS: "<<M_OS<<", M_TM: "<<M_TM<<endl;
	cout<<"SL/M_OS: "<<djack_t(SL_OS/M_OS)<<", SL/M_TM: "<<djack_t(SL_TM/M_TM)<<endl;
	cout<<"DP_OS: "<<DZP_fr_ZP_OS<<", DP_TM: "<<DZP_fr_ZP_TM<<endl;
	cout<<"DA0_OS: "<<DZA_fr_ZA_OS<<", DA0_TM: "<<DZA_fr_ZA_TM<<endl;

	}
	djack_t za_QED_fact_exp=DZP_fr_ZP_TM+SL_TM/M_TM-DZA_fr_ZA_OS+use_mass*pred_Zp_coeff/sqr(4*M_PI);
	djack_t zv_QED_fact_exp=DZP_fr_ZP_TM+SL_TM/M_TM-DZA_fr_ZA_TM+use_mass*pred_Zp_coeff/sqr(4*M_PI);
	djack_t za_m_zv_QED_fact_exp=DZA_fr_ZA_TM-DZA_fr_ZA_OS;
	djack_t zp_m_zs_QED_fact_exp=DZP_fr_ZP_OS-DZP_fr_ZP_TM;
	
	//cout<<"Factorization Za expanded from dec: "<<smart_print(djack_t(za_QED_fact_exp_dec*e2).ave_err())<<endl;
	// cout<<"Factorization Za expanded: "<<smart_print(djack_t(za_QED_fact_exp*e2).ave_err())<<endl;
	// cout<<"Factorization Zv expanded: "<<smart_print(djack_t(zv_QED_fact_exp*e2).ave_err())<<endl;
	djack_t zp_fr_zs=ZP_OS/ZP_TM;
	djack_t zv_fr_za=ZA_TM/ZA_OS;
	djack_t za_QED_fact=za_QED_fact_exp*sqr(4*M_PI);
	djack_t zv_QED_fact=zv_QED_fact_exp*sqr(4*M_PI);
	djack_t zp_m_zs_QED_fact=zp_m_zs_QED_fact_exp*sqr(4*M_PI);
	djack_t za_m_zv_QED_fact=za_m_zv_QED_fact_exp*sqr(4*M_PI);
	djack_t za_QED_fact_corr=za_QED_fact/pred_Za_coeff-1.0;
	djack_t zv_QED_fact_corr=zv_QED_fact/pred_Zv_coeff-1.0;
	djack_t zp_m_zs_QED_fact_corr=zp_m_zs_QED_fact/pred_Zp_m_Zs_coeff-1.0;
	djack_t za_m_zv_QED_fact_corr=za_m_zv_QED_fact/pred_Za_m_Zv_coeff-1.0;
	cout<<"Zp/Zs: "<<smart_print(zp_fr_zs)<<endl;
	cout<<"Zv/Za: "<<smart_print(zv_fr_za)<<endl;
	cout<<"Factorization Za: "<<smart_print(za_QED_fact)<<", Aoki prediction: "<<pred_Za_coeff<<", correction: "<<smart_print(za_QED_fact_corr)<<endl;
	cout<<"Factorization Zv: "<<smart_print(zv_QED_fact)<<", Aoki prediction: "<<pred_Zv_coeff<<", correction: "<<smart_print(zv_QED_fact_corr)<<endl;
	cout<<"Factorization Zp-Zs: "<<smart_print(zp_m_zs_QED_fact)<<", Aoki prediction: "<<pred_Zp_m_Zs_coeff<<", correction: "<<smart_print(zp_m_zs_QED_fact_corr)<<endl;
	cout<<"Factorization Za-Zv: "<<smart_print(za_m_zv_QED_fact)<<", Aoki prediction: "<<pred_Za_m_Zv_coeff<<", correction: "<<smart_print(za_m_zv_QED_fact_corr)<<endl;
      }
      
      //determination of Za-Zv from vector corr
      bool VEC=false;
      if(VEC)
      {
	cout<<endl<<"Determining Za-ZV from vector correlator"<<endl;
	
	auto load_LO_CORR=[&e2_var,&ka_var,&ma_var,&e2_phys,&deltam_cr,&deltam,&use_mass,&ant]
	  (djvec_t &LO,djvec_t &DELTA,djvec_t &CORR,const string &name,const int &rcombo,const int &reim,const int &tpar,const int &rpar)
	  {
	    LO=load_2pts_spave(name,rcombo,reim,tpar,rpar);
	    LO.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_LO.xmg");
	    djvec_t E2;
	    if(ant) E2=der_2pts_spave(name,rcombo,reim,tpar,rpar,e2_var,"E2E2");
	    else    E2=
		      load_2pts_spave(name+"_LL",rcombo,reim,tpar,rpar)+
		      load_2pts_spave(name+"_0M",rcombo,reim,tpar,rpar)+
		      load_2pts_spave(name+"_M0",rcombo,reim,tpar,rpar)+
		      load_2pts_spave(name+"_0T",rcombo,reim,tpar,rpar)+
		      load_2pts_spave(name+"_T0",rcombo,reim,tpar,rpar);
	    E2.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_E2.xmg");
	    
	    djvec_t KA;
	    if(ant) KA=der_2pts_spave(name,rcombo,reim,tpar,rpar,ka_var,"KAKA");
	    else    KA=-2*load_2pts_spave(name+"_P0",rcombo,!reim,tpar,!rpar);
	    KA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_KA.xmg");
	    
	    djvec_t MA;
	    if(ant) MA=der_2pts_spave(name,rcombo,reim,tpar,rpar,ma_var,"MAMA");
	    else    MA=-(load_2pts_spave(name+"_0S",rcombo,reim,tpar,rpar)+load_2pts_spave(name+"_S0",rcombo,reim,tpar,rpar));
	    MA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_MA.xmg");
	    
	    DELTA=E2-deltam_cr*KA+deltam*MA*use_mass;
	    DELTA.ave_err().write("plots/"+name+"_"+RTAGS_2PTS_NAME[rcombo]+"_DELTA.xmg");
	    CORR=LO+e2_phys*DELTA;
	  };
	
	djvec_t LO_VKVK_SAME,DELTA_VKVK_SAME,CORR_VKVK_SAME;
	load_LO_CORR(LO_VKVK_SAME,DELTA_VKVK_SAME,CORR_VKVK_SAME,"V%cV%c",RSAME,RE,EVN,EVN);
	djvec_t LO_VKVK_OPPO,DELTA_VKVK_OPPO,CORR_VKVK_OPPO;
	load_LO_CORR(LO_VKVK_OPPO,DELTA_VKVK_OPPO,CORR_VKVK_OPPO,"V%cV%c",ROPPO,RE,EVN,EVN);
	
	djack_t M_OS,SL_OS;
	djack_t Z2_OS,DZ2_fr_Z2_OS;
	two_pts_with_ins_ratio_fit(Z2_OS,M_OS,DZ2_fr_Z2_OS,SL_OS,LO_VKVK_OPPO,DELTA_VKVK_OPPO,T/2,tmin,tmax,"plots/VKVK_OPPO_effmass.xmg","plots/VKVK_OPPO_slope.xmg");
	djack_t ZV_OS=sqrt(Z2_OS);
	djack_t DZV_fr_ZV_OS=DZ2_fr_Z2_OS/2;
	
	djack_t M_TM,SL_TM;
	djack_t Z2_TM,DZ2_fr_Z2_TM;
	two_pts_with_ins_ratio_fit(Z2_TM,M_TM,DZ2_fr_Z2_TM,SL_TM,LO_VKVK_SAME,DELTA_VKVK_SAME,T/2,tmin,tmax,"plots/VKVK_SAME_effmass.xmg","plots/VKVK_SAME_slope.xmg");
	djack_t ZV_TM=sqrt(Z2_TM);
	djack_t DZV_fr_ZV_TM=DZ2_fr_Z2_TM/2;
	
	cout<<"SL_OS: "<<SL_OS<<", SL_TM: "<<SL_TM<<endl;
	cout<<"M_OS: "<<M_OS<<", M_TM: "<<M_TM<<endl;
	cout<<"SL/M_OS: "<<djack_t(SL_OS/M_OS)<<", SL/M_TM: "<<djack_t(SL_TM/M_TM)<<endl;
	cout<<"DV_OS: "<<DZV_fr_ZV_OS<<", DV_TM: "<<DZV_fr_ZV_TM<<endl;
	
	djack_t za_m_zv_QED_fact_exp=DZV_fr_ZV_OS-DZV_fr_ZV_TM;
	djack_t zv_fr_za=ZV_OS/ZV_TM;
	djack_t za_m_zv_QED_fact=za_m_zv_QED_fact_exp*sqr(4*M_PI);
	djack_t za_m_zv_QED_fact_corr=za_m_zv_QED_fact/pred_Za_m_Zv_coeff-1.0;
	cout<<"Zv/Za: "<<smart_print(zv_fr_za)<<endl;
	cout<<"Factorization Za-Zv: "<<smart_print(za_m_zv_QED_fact)<<", Aoki prediction: "<<pred_Za_m_Zv_coeff<<", correction: "<<smart_print(za_m_zv_QED_fact_corr)<<endl;
      }
      
      bool WI_cons_cur=false;
      if(WI_cons_cur)
	{
	  //zv non conserved
	  djvec_t LO_VKVK_zv=load_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN);
	  djvec_t E2_VKVK_zv=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,e2_var,"E2E2");
	  djvec_t KA_VKVK_zv=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,ka_var,"KAKA");
	  djvec_t CORR_VKVK_zv=LO_VKVK_zv+e2_phys*djvec_t(E2_VKVK_zv-deltam_cr*KA_VKVK_zv);
	  djack_t Z2_zv,M_zv;
	  two_pts_fit(Z2_zv,M_zv,LO_VKVK_zv,T/2,tmin,tmax,"plots/effmass_M_VKVK_zv.xmg");
	  djack_t CORR_Z2_zv,CORR_M_zv;
	  two_pts_fit(CORR_Z2_zv,CORR_M_zv,CORR_VKVK_zv,T/2,tmin,tmax,"plots/effmass_M_CORR_VKVK_zv.xmg");
	  
	  //za non conserved, to be matched
	  djvec_t LO_VKVK_za=load_2pts_spave("V%cV%c",RSAME,RE,EVN,EVN);
	  djvec_t E2_VKVK_za=der_2pts_spave("V%cV%c",RSAME,RE,EVN,EVN,e2_var,"E2E2");
	  djvec_t KA_VKVK_za=der_2pts_spave("V%cV%c",RSAME,RE,EVN,EVN,ka_var,"KAKA");
	  djvec_t CORR_VKVK_za=LO_VKVK_za+e2_phys*djvec_t(E2_VKVK_za-deltam_cr*KA_VKVK_za);
	  djack_t Z2_za,M_za;
	  two_pts_fit(Z2_za,M_za,LO_VKVK_za,T/2,tmin,tmax,"plots/effmass_M_VKVK_za.xmg");
	  djack_t CORR_Z2_za,CORR_M_za;
	  two_pts_fit(CORR_Z2_za,CORR_M_za,CORR_VKVK_za,T/2,tmin,tmax,"plots/effmass_M_CORR_VKVK_za.xmg");
	  
	  //compute zv/za from direct ratio
	  grace_file_t zv_fr_za_t_file("plots/LO_VKVK_zv_fr_za.xmg");
	  djvec_t zv_fr_za_t=djvec_t(sqrt(LO_VKVK_za/LO_VKVK_zv));
	  zv_fr_za_t_file.write_vec_ave_err(zv_fr_za_t.ave_err());
	  djack_t zv_fr_za=constant_fit(zv_fr_za_t,tmin,tmax);
	  zv_fr_za_t_file.write_constant_band(tmin,tmax,zv_fr_za);
	  djack_t zv_fr_za_ri=djack_t(gauss_filler_t{Zv_ae[0][ib],1234})/djack_t(gauss_filler_t{Za_ae[0][ib],1234});
	  zv_fr_za_t_file.write_constant_band(tmin,tmax,zv_fr_za_ri);
	  
	  //compute corr to zv/za from direct ratio
	  grace_file_t CORR_zv_fr_za_t_file("plots/CORR_VKVK_zv_fr_za.xmg");
	  djvec_t CORR_zv_fr_za_t=djvec_t(sqrt(CORR_VKVK_za/CORR_VKVK_zv));
	  CORR_zv_fr_za_t_file.write_vec_ave_err(CORR_zv_fr_za_t.ave_err());
	  djack_t CORR_zv_fr_za=constant_fit(CORR_zv_fr_za_t,tmin,tmax);
	  CORR_zv_fr_za_t_file.write_constant_band(tmin,tmax,CORR_zv_fr_za);
	  
	  cout<<"Zv/Za:      "<<smart_print(djack_t(sqrt(Z2_za/Z2_zv)).ave_err())<<endl;
	  cout<<"CORR Zv/Za: "<<smart_print(djack_t(sqrt(CORR_Z2_za/CORR_Z2_zv)).ave_err())<<endl;
	  cout<<"DELTA Zv-Za: "<<smart_print(djack_t(sqrt(CORR_Z2_za/CORR_Z2_zv)-sqrt(Z2_za/Z2_zv)).ave_err())<<endl;
	  
	  //zv conserved
	  djvec_t LO_VKCVK_zv=load_2pts_spave("S0V%c",ROPPO,RE,EVN,EVN,"","","V%c");
	  djvec_t E2_VKCVK_zv=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,e2_var,"E2","","V%c");
	  djvec_t KA_VKCVK_zv=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,ka_var,"KA","","V%c");
	  djvec_t CORR_VKCVK_zv=LO_VKCVK_zv+e2_phys*djvec_t(E2_VKCVK_zv-deltam_cr*KA_VKCVK_zv);
	  
	  djvec_t ZV_ren_t=LO_VKCVK_zv/LO_VKVK_zv;
	  djvec_t CORR_ZV_ren_t=CORR_VKCVK_zv/CORR_VKVK_zv;
	  djack_t ZV_ren=constant_fit(ZV_ren_t,tmin,tmax,"plots/ZV_WI_2pts.xmg");
	  djack_t CORR_ZV_ren=constant_fit(CORR_ZV_ren_t,tmin,tmax,"plots/CORR_ZV_WI_2pts.xmg");
	  djvec_t ZV_QED_fact_t=(CORR_ZV_ren_t-ZV_ren_t)/(ZV_ren_t*qf2_phys);
	  djack_t ZV_QED_fact=constant_fit(ZV_QED_fact_t,T/4-2,T/4+2,"plots/ZV_WI_2pts_QED_fact.xmg");
	  
	  djvec_t LO_VKVK=load_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN);
	  djvec_t E2_VKVK=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,e2_var,"E2");
	  djvec_t KA_VKVK=der_2pts_spave("V%cV%c",ROPPO,RE,EVN,EVN,ka_var,"KA");
	  djvec_t CORR_VKVK_v=LO_VKVK+e2_phys*djvec_t(E2_VKVK-deltam_cr*KA_VKVK);
	  
	  //compute fact Za
	  djvec_t ZA_ren_t=ZV_ren_t/zv_fr_za_t;
	  djvec_t CORR_ZA_ren_t=CORR_ZV_ren_t/CORR_zv_fr_za_t;
	  djvec_t ZA_QED_fact_t=(CORR_ZA_ren_t-ZA_ren_t)/(ZA_ren_t*qf2_phys);
	  djack_t ZA_QED_fact=constant_fit(ZA_QED_fact_t,T/4-2,T/4+2,"plots/ZA_WI_2pts_QED_fact.xmg");
	}
    }
  
  return 0;
}
