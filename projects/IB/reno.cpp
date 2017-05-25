#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <common.hpp>

const int EVN=1,ODD=-1;//,UNK=0;

typedef int rcombo_2pts_t[2][2];
typedef int rcombo_3pts_t[2][3];

const int SAMER_2PTS=0;//,OPPOR_2PTS=1;
const rcombo_2pts_t RTAGS_2PTS[2]={{{0,0},{1,1}}, {{0,1},{1,0}}};
const array<string,2> RTAGS_2PTS_NAME={"SAME","OPPO"};

//const int SAMER_SAMER_3PTS=0,SAMER_OPPOR_3PTS=1,OPPOR_SAMER_3PTS=2,OPPOR_OPPOR_3PTS=3;
const array<string,4> RTAGS_3PTS_NAME={"SAME_SAME","SAME_OPPO","OPPO_SAME","OPPO_OPPO"};
const rcombo_3pts_t RTAGS_3PTS[4]={{{0,0,1},{1,1,0}}, {{0,0,0},{1,1,1}}, {{1,0,1},{0,1,0}}, {{0,1,1},{1,0,0}}}; //REV, SPEC, SEQ , remember that SEQ is already reversed

size_t T;
double amq;
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
{return Zso*Zsi/(4*Mso*Msi)*exp(-Mso*t)*exp(-Msi*(T/2-t))*(Mso+Msi);}

//! compute the Z from 2pts decay constant
void z_from_2pts_dec(djack_t &za,djack_t &zv,const djvec_t &P5P5_SAME,const djvec_t &A0P5_SAME,const djvec_t &P5P5_OPPO,const djvec_t &A0P5_OPPO,const string &tag)
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
  djack_t f_SAME=ZP_SAME*(2*amq)/sqr(MP5P5_SAME);
  djack_t f_SAME_bare=ZA_SAME/MP5P5_SAME;
  
  djack_t ZA_OPPO=ZA0P5_OPPO/sqrt(ZP5P5_OPPO);
  djack_t f_OPPO_bare=ZA_OPPO/MP5P5_OPPO;
  
  za=f_SAME/f_OPPO_bare;
  zv=f_SAME/f_SAME_bare;
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
  
  amq=input.read<double>("amq");
  const double k=input.read<double>("k");
  const double kp=input.read<double>("kp");
  const double km=input.read<double>("km");
  const size_t ib=input.read<size_t>("ib");
  
  const size_t WI_2pts=input.read<size_t>("WI_2pts");
  
  const double qf2_der=sqr(0.25/3),qf2_phys=1.0/9.0;
  const array<double,3> e2_var={0,qf2_der*e2,qf2_der*e2};
  const double e2_phys=e2*qf2_phys;
  
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  const array<double,3> ka_var={mc,mcp,mcm};
  
  //compute deltamcr
  djvec_t E2_V0P5=der_2pts("V0P5",SAMER_2PTS,IM,ODD,ODD,e2_var,"E2");
  djvec_t KA_V0P5=der_2pts("V0P5",SAMER_2PTS,IM,ODD,ODD,ka_var,"KA");
  djvec_t deltam_cr_t=E2_V0P5/KA_V0P5;
  deltam_cr_t.ave_err().write("plots/deltam_cr.xmg");
  djack_t deltam_cr=deltam_cr_t[13];
  
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
  
  const double pred_Zv=-20.6178*alpha_em/(4*M_PI);
  const double pred_Za=-15.7963*alpha_em/(4*M_PI);
  
  for(int rdiff_so=0;rdiff_so<2;rdiff_so++)
    for(int rdiff_si=0;rdiff_si<2;rdiff_si++)
      {
	int rdiff_tot=rdiff_si+2*rdiff_so;
	cout<<"Case: "<<RTAGS_3PTS_NAME[rdiff_tot]<<endl;
	
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
		(LO_P5P5[rdiff_so][t]-two_pts_corr_fun(Z2[rdiff_so],M[rdiff_so],T/2,T-t,0))*
		(LO_P5P5[rdiff_si][T/2-t]-two_pts_corr_fun(Z2[rdiff_si],M[rdiff_si],T/2,T/2+t,0))
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
	    eff_file.write_constant_band(0,T/2,djack_t(M[rdiff_so]-M[rdiff_si]));
	  }
	  
	  //write the aperiodic effective mass of 3pts corrected minus non corrected
	  {
	    grace_file_t eff_corr_file("plots/effmass_3pts_corr_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	    eff_corr_file.write_vec_ave_err(djvec_t(effective_mass(CORR_3,T/2,0,1)-effective_mass(LO_3,T/2,0,1)).ave_err());
	    eff_corr_file.write_constant_band(0,T/2,djack_t(M_CORR[rdiff_so]-M_CORR[rdiff_si]-M[rdiff_so]+M[rdiff_si]));
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
  
  if(WI_2pts)
    {
      int RSAME=0,ROPPO=1;
      
      //determination of Za from 2pts WI
      {
	cout<<"Determining Za and Zav from 2pts decay constant"<<endl;
	
	auto load_LO_CORR=[&e2_var,&ka_var,&e2_phys,&deltam_cr]
	  (djvec_t &LO,djvec_t &CORR,const string &name,const int &rcombo,const int &reim,const int &tpar,const int &rpar)
	  {
	    LO=load_2pts(name,rcombo,reim,tpar,rpar);
	    djvec_t E2=der_2pts(name,rcombo,reim,tpar,rpar,e2_var,"E2E2");
	    djvec_t KA=der_2pts(name,rcombo,reim,tpar,rpar,ka_var,"KAKA");
	    CORR=LO+e2_phys*djvec_t(E2-deltam_cr*KA);
	  };
	
	djvec_t LO_P5P5_SAME,CORR_P5P5_SAME;
	load_LO_CORR(LO_P5P5_SAME,CORR_P5P5_SAME,"P5P5",RSAME,RE,EVN,EVN);
	djvec_t LO_P5P5_OPPO,CORR_P5P5_OPPO;
	load_LO_CORR(LO_P5P5_OPPO,CORR_P5P5_OPPO,"P5P5",ROPPO,RE,EVN,EVN);
	djvec_t LO_A0P5_SAME,CORR_A0P5_SAME;
	load_LO_CORR(LO_A0P5_SAME,CORR_A0P5_SAME,"A0P5",RSAME,RE,ODD,EVN);
	djvec_t LO_A0P5_OPPO,CORR_A0P5_OPPO;
	load_LO_CORR(LO_A0P5_OPPO,CORR_A0P5_OPPO,"A0P5",ROPPO,RE,ODD,EVN);
	
	djack_t za,zv;
	z_from_2pts_dec(za,zv,LO_P5P5_SAME,LO_A0P5_SAME,LO_P5P5_OPPO,LO_A0P5_OPPO,"LO");
	
	cout<<"Za: "<<za.ave_err()<<endl;
	cout<<"Zv: "<<zv.ave_err()<<endl;
	
	djack_t CORR_za,CORR_zv;
	z_from_2pts_dec(CORR_za,CORR_zv,CORR_P5P5_SAME,CORR_A0P5_SAME,CORR_P5P5_OPPO,CORR_A0P5_OPPO,"CORR");
	
	cout<<"Corr Za: "<<CORR_za.ave_err()<<endl;
	cout<<"Corr Zv: "<<CORR_zv.ave_err()<<endl;
	
	djack_t zv_QED_fact=(CORR_zv-zv)/(zv*qf2_phys);
	djack_t za_QED_fact=(CORR_za-za)/(za*qf2_phys);
      
	cout<<"Factorization Za: "<<za_QED_fact.ave_err()<<endl;
	cout<<"Factorization Zv: "<<zv_QED_fact.ave_err()<<endl;
      }
      
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
      djack_t zv_fr_za_ri=djack_t({Zv_ae[0][ib],1234})/djack_t({Za_ae[0][ib],1234});
      zv_fr_za_t_file.write_constant_band(tmin,tmax,zv_fr_za_ri);
      
      //compute corr to zv/za from direct ratio
      grace_file_t CORR_zv_fr_za_t_file("plots/CORR_VKVK_zv_fr_za.xmg");
      djvec_t CORR_zv_fr_za_t=djvec_t(sqrt(CORR_VKVK_za/CORR_VKVK_zv));
      CORR_zv_fr_za_t_file.write_vec_ave_err(CORR_zv_fr_za_t.ave_err());
      djack_t CORR_zv_fr_za=constant_fit(CORR_zv_fr_za_t,tmin,tmax);
      CORR_zv_fr_za_t_file.write_constant_band(tmin,tmax,CORR_zv_fr_za);
      
      cout<<"Zv/Za:      "<<smart_print(djack_t(sqrt(Z2_za/Z2_zv)).ave_err())<<endl;
      cout<<"CORR Zv/Za: "<<smart_print(djack_t(sqrt(CORR_Z2_za/CORR_Z2_zv)).ave_err())<<endl;
      
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
  
  return 0;
}
