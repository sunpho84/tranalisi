#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const int EVN=1,UNK=0,ODD=-1;

typedef int rcombo_2pts_t[2][2];
typedef int rcombo_3pts_t[2][3];

const int SAMER_2PTS=0,OPPOR_2PTS=1;
const rcombo_2pts_t RTAGS_2PTS[2]={{{0,0},{1,1}}, {{0,1},{1,0}}};
const array<string,2> RTAGS_2PTS_NAME={"SAME","OPPO"};

const int SAMER_SAMER_3PTS=0,SAMER_OPPOR_3PTS=1,OPPOR_SAMER_3PTS=2,OPPOR_OPPOR_3PTS=3;
const array<string,4> RTAGS_3PTS_NAME={"SAME_SAME","SAME_OPPO","OPPO_SAME","OPPO_OPPO"};
const rcombo_3pts_t RTAGS_3PTS[4]={{{0,0,1},{1,1,0}}, {{0,0,0},{1,1,1}}, {{1,0,1},{0,1,0}}, {{0,1,1},{1,0,0}}}; //REV, SPEC, SEQ , remember that SEQ is already reversed

const size_t T=48;
  const int tmin=12,tmax=T/2-2;

//! return "_var" or ""
string prespaced(const string &var)
{return (var==""?"":"_")+var;}

//! read a single correlator
djvec_t get(const string &name,const int ri,const int tpar)
{
  cout<<"Reading "<<name<<endl;
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
djvec_t load_averaging(const array<string,2> &what,int ri,int tpar,int rpar,const string &diracs,const string &var="")
{
  djvec_t temp[2];
  for(int r=0;r<2;r++)
    temp[r]=get(what[r]+prespaced(var)+"_"+diracs,ri,tpar);
  
  return djvec_t(temp[0]+rpar*temp[1])/(1+abs(rpar));
}

//! load 2pts (??P5)
djvec_t load_2pts(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const string &var="")
{
  array<string,2> what;
  const rcombo_2pts_t &rtags=RTAGS_2PTS[irtags];
  for(int r=0;r<2;r++) what[r]="Spect"+to_string(rtags[r][0])+"_Spect"+to_string(rtags[r][1]);
  
  djvec_t out=load_averaging(what,ri,tpar,rpar,dirtag+"P5",var);
  out.ave_err().write("plots/"+dirtag+"P5_"+RTAGS_2PTS_NAME[irtags]+prespaced(var)+".xmg");
  
  return out;
}

//! compute the derivative od 2pts correlators w.r.t a given variation
djvec_t der_2pts(const string &dirtag,const int irtags,const int ri,const int tpar,const int rpar,const array<double,3> &coeffs,const string &var)
{
  array<string,3> vars={"",var+"P",var+"M"};
  array<djvec_t,3> corrs;
  for(int i=0;i<3;i++) corrs[i]=load_2pts(dirtag,irtags,ri,tpar,rpar,vars[i]);
  
  djvec_t out=der(corrs,coeffs);
  out.ave_err().write("plots/"+dirtag+"P5_"+RTAGS_2PTS_NAME[irtags]+"_der_"+var+".xmg");
  
  return out;
}

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

int main()
{
  set_njacks(15);
  
  const double qf2_der=sqr(0.25/3),qf2_phys=1.0/9.0;
  const array<double,3> e2_var={0,qf2_der*e2,qf2_der*e2};
  const double e2_phys=e2*qf2_phys;
  
  const double k=0.163255,kp=0.163305,km=0.163205;
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  const array<double,3> ka_var={mc,mcp,mcm};
  
  //compute deltamcr
  djvec_t E2_V0P5=der_2pts("V0",SAMER_2PTS,IM,ODD,ODD,e2_var,"E2");
  djvec_t KA_V0P5=der_2pts("V0",SAMER_2PTS,IM,ODD,ODD,ka_var,"KA");
  djvec_t deltam_cr_corr=E2_V0P5/KA_V0P5;
  deltam_cr_corr.ave_err().write("plots/deltam_cr.xmg");
  djack_t deltam_cr=deltam_cr_corr[tmin];
  
  djvec_t Z(2),Z2(2),M(2);
  djvec_t Z_CORR(2),Z2_CORR(2),M_CORR(2);
  djvec_t LO_P5P5[2];
  for(int rdiff=0;rdiff<2;rdiff++)
    {
      LO_P5P5[rdiff]=load_2pts("P5",rdiff,RE,EVN,EVN);
      two_pts_fit(Z2[rdiff],M[rdiff],LO_P5P5[rdiff],T/2,tmin,tmax,"plots/effmass_P5P5_"+RTAGS_2PTS_NAME[rdiff]+".xmg");
      Z[rdiff]=sqrt(Z2[rdiff]);
      
      djvec_t E2_P5P5=der_2pts("P5",rdiff,RE,EVN,EVN,e2_var,"E2");
      djvec_t KA_P5P5=der_2pts("P5",rdiff,RE,EVN,EVN,ka_var,"KA");
      
      djvec_t CORR_P5P5=LO_P5P5[rdiff]+e2_phys*djvec_t(E2_P5P5-deltam_cr*KA_P5P5);
      CORR_P5P5.ave_err().write("plots/P5P5_"+RTAGS_2PTS_NAME[rdiff]+"_corrected.xmg");
      two_pts_fit(Z2_CORR[rdiff],M_CORR[rdiff],CORR_P5P5,T/2,tmin,tmax,"plots/effmass_P5P5_"+RTAGS_2PTS_NAME[rdiff]+"_corrected.xmg");
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
	
	djvec_t LO_3=load_3pts(rdiff_tot,ODD,EVN);
	
	//LO
	djack_t Zren;
	{
	  //reconstruct 2pts
	  djvec_t LO_2_reco(T/2);
	  for(size_t t=0;t<T/2;t++) LO_2_reco[t]=Z[rdiff_so]*Z[rdiff_si]/(4*M[rdiff_so]*M[rdiff_si])*exp(-M[rdiff_so]*t)*exp(-M[rdiff_si]*(T/2-t));
	  cout<<"Check "<<RTAGS_3PTS_NAME[rdiff_tot]<<"reco, "<<LO_2_reco[0].ave_err()<<" "<<djack_t(LO_P5P5[rdiff_so][T/2]/djack_t(2*M[rdiff_so]+2*M[rdiff_si])).ave_err()<<endl;
	  
	  djvec_t INVMATREL=-LO_2_reco/LO_3,Zren_corr=INVMATREL*djack_t(M[rdiff_so]+M[rdiff_si]);
	  Zren_corr.ave_err().write("plots/Z_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	  Zren=Zren_corr[T/4];
	}
	
	//QED
	djack_t CORR_Zren;
	{
	  //build corrected 3pts
	  djvec_t E2_3=der_3pts(rdiff_tot,ODD,EVN,e2_var,"E2");
	  djvec_t KA_3=der_3pts(rdiff_tot,ODD,EVN,ka_var,"KA");
	  djvec_t CORR_3=LO_3+e2_phys*djvec_t(E2_3-deltam_cr*KA_3);
	  
	  //reconstruct 2pts corrected
	  djvec_t CORR_2_reco(T/2);
	  for(size_t t=0;t<T/2;t++) CORR_2_reco[t]=Z_CORR[rdiff_so]*Z_CORR[rdiff_si]/(4*M_CORR[rdiff_so]*M_CORR[rdiff_si])*exp(-M_CORR[rdiff_so]*t)*exp(-M_CORR[rdiff_si]*(T/2-t));
	  CORR_2_reco.ave_err().write("plots/CORR_2_reco_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	  djvec_t CORR_INVMATREL=-CORR_2_reco/CORR_3,CORR_Zren_corr=CORR_INVMATREL*djack_t(M_CORR[rdiff_so]+M_CORR[rdiff_si]);
	  CORR_Zren_corr.ave_err().write("plots/CORR_Z_"+RTAGS_3PTS_NAME[rdiff_tot]+".xmg");
	  CORR_Zren=CORR_Zren_corr[T/4];
	}
	
	djack_t fact=(CORR_Zren-Zren)/(Zren*qf2_phys);
	cout<<"Factorization: "<<fact.ave_err()<<endl;
      }
  
  return 0;
}
