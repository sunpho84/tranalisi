#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=48;
const size_t nconfs=150;
const range_t conf_set={700,10,700+10*(nconfs-1)};

//! return the whole correlator
djvec_t get(const string &combo,const string &ID,const size_t reim,const size_t sym)
{
  const vector<string> cID=
    {"S0P5","V1P5","V2P5","V3P5","V0P5","P5P5","A1P5","A2P5","A3P5","A0P5","T1P5","T2P5","T3P5","B1P5","B2P5","B3P5","P5S0","P5V1","P5V2","P5V3","P5V0",
     "P5P5","P5A1","P5A2","P5A3","P5A0","P5T1","P5T2","P5T3","P5B1","P5B2","P5B3","P5P5"};
 
 const size_t ncols=2;
 const string templ="data/%04d/mes_contr_"+combo;
 const djvec_t data=read_conf_set_t(templ,conf_set,ncols,{0,1},T,SILENT);
 if(data.size()==0) CRASH("No file opened for template %s",templ.c_str());
 
 const auto pos=find(cID.begin(),cID.end(),ID);
 if(pos==cID.end()) CRASH("Unknown %s",ID.c_str());
 
 //filter
 const size_t offset=ncols*distance(cID.begin(),pos)+reim;
 const size_t each=ncols*cID.size();
 const size_t base_nel=T;
 const size_t hw=data.size()/(base_nel*each);
 djvec_t res=vec_filter(data,gslice(base_nel*offset,{hw,T},{each*base_nel,1}));
 
 //write
 res.ave_err().write("plots/"+combo+"_"+ID+".xmg");
 
 //self-symmetrize
 res.symmetrize(sym);
 
 //write symmetrized
 res.ave_err().write("plots/"+combo+"_"+ID+"_symm.xmg");
 
 return res;
}

djvec_t compute_der(const string &qrev,const string &qfo,const string &PV_rev,const string &PV_fo
		    ,const string &MV_rev,const string &MV_fo,const string &ID,const size_t reim,const size_t sym,const double varP,const double varM,const string &name)
{
  djvec_t X=get(qrev+"_"+qfo,ID,reim,sym);
  djvec_t XP=get(qrev+PV_rev+"_"+qfo+PV_fo,ID,reim,sym);
  djvec_t XM=get(qrev+MV_rev+"_"+qfo+MV_fo,ID,reim,sym);
  djvec_t X_der_P=djvec_t(XP-X)/varP;
  X_der_P.ave_err().write("plots/"+qrev+"_"+qfo+"_"+ID+"_der_"+name+"P.xmg");
  djvec_t X_der_M=djvec_t(XM-X)/varM;
  X_der_M.ave_err().write("plots/"+qrev+"_"+qfo+"_"+ID+"_der_"+name+"M.xmg");
  djvec_t X_der=djvec_t(X_der_P+X_der_M)/2;
  X_der.ave_err().write("plots/"+qrev+"_"+qfo+"_"+ID+"_der_"+name+".xmg");
  
  return X_der;
}

int main()
{
  const double m=0.02363,mp=0.02373,mm=0.02353;
  const double k=0.163255,kp=0.163355,km=0.163155;
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  const double qf2_der=sqr(0.25/3);
  
  set_njacks(15);
  
  //2pts
  djvec_t Spect_Spect=get("Spect_Spect","P5P5",RE,1);
  djack_t M,Z2;
  two_pts_fit(Z2,M,Spect_Spect,T/2,12,T/2,"plots/Spect_Spect_effmass.xmg");
  cout<<"Mass: "<<M.ave_err()<<endl;
  
  djvec_t Spect_Spect_der_e2=compute_der("Spect","Spect","P","","M","","P5P5",RE,1,qf2_der*e2,qf2_der*e2,"e2");
  djvec_t deltam_cr_num=compute_der("Spect","Spect","P","","M","","V0P5",IM,-1,qf2_der*e2,qf2_der*e2,"e2");
  djvec_t deltam_cr_den=compute_der("Spect","Spect","KAP","","KAM","","V0P5",IM,-1,mcp-mc,mcm-mc,"KA");
  djvec_t deltam_cr_corr=deltam_cr_num/deltam_cr_den;
  deltam_cr_corr.ave_err().write("plots/deltam_cr.xmg");
  djack_t deltam_cr=deltam_cr_corr[T/4];
  djvec_t Spect_Spect_der_MA=compute_der("Spect","Spect","MAP","","MAM","","P5P5",RE,1,mp-m,mm-m,"MA");
  djvec_t Spect_Spect_der_KA=compute_der("Spect","Spect","KAP","","KAM","","P5P5",RE,1,mcp-mc,mcm-mc,"KA");
  
  //3pts
  djvec_t Spect_Seq=get("Spect_Seq","V0P5",RE,-1);
  djvec_t Spect_V_Seq=get("Spect_V_Seq","S0P5",RE,-1);
  
  djvec_t Spect_Seq_der_e2=compute_der("Spect","Seq","P","P","M","M","V0P5",RE,-1,qf2_der*e2,qf2_der*e2,"e2");
  djvec_t Spect_Seq_der_MA=compute_der("Spect","Seq","MAP","MAP","MAM","MAM","V0P5",RE,-1,mp-m,mm-m,"MA");
  djvec_t Spect_Seq_der_KA=compute_der("Spect","Seq","KAP","KAP","KAM","KAM","V0P5",RE,-1,mcp-mc,mcm-mc,"KA");
  
  djvec_t Spect_V_Seq_der_e2=compute_der("Spect","V_Seq","P","P","M","M","S0P5",RE,-1,qf2_der*e2,qf2_der*e2,"e2");
  djvec_t Spect_V_Seq_der_MA=compute_der("Spect","V_Seq","MAP","MAP","MAM","MAM","S0P5",RE,-1,mp-m,mm-m,"MA");
  djvec_t Spect_V_Seq_der_KA=compute_der("Spect","V_Seq","KAP","KAP","KAM","KAM","S0P5",RE,-1,mcp-mc,mcm-mc,"KA");
  
  //Zv
  cout<<"Zv_num: "<<Spect_Spect[T/2].ave_err()<<endl;
  djvec_t Zv=Spect_Spect[T/2]/(-2.0*Spect_Seq);
  Zv.ave_err().write("plots/Zv.xmg");
  djvec_t CZv=Spect_Spect[T/2]/(-2.0*Spect_V_Seq);
  CZv.ave_err().write("plots/CZv.xmg");
  
  //3pts derivative w.r.t to e2
  // djvec_t SpectP_CV_SeqP=get("SpectP_V_SeqP","V0P5",RE,1);
  // djvec_t SpectM_CV_SeqM=get("SpectM_V_SeqM","V0P5",RE,1);
  // djvec_t Spect_CV_Seq_der_P=djvec_t(SpectP_CV_SeqP-Spect_CV_Seq)/(qf2_der*e2);
  // Spect_CV_Seq_der_P.ave_err().write("plots/Spect_CV_Seq_der_P.xmg");
  // djvec_t Spect_CV_Seq_der_M=djvec_t(SpectM_CV_SeqM-Spect_CV_Seq)/(qf2_der*e2);
  // Spect_CV_Seq_der_M.ave_err().write("plots/Spect_CV_Seq_der_M.xmg");
  // djvec_t Spect_CV_Seq_der_=djvec_t(Spect_CV_Seq_der_P+Spect_CV_Seq_der_M)/2;
  // Spect_CV_Seq_der_.ave_err().write("plots/Spect_CV_Seq_der_.xmg");
  
  //QED_corr to Zv
  double amq=0.02363;
  double a=0.4497;
  double deltam=amq*(6.0*log(mu_MS*a)-22.596)/(16.0*sqr(M_PI));
  
  //deltam is uselsess, begins at higher orders
  deltam*=1.0;
  
  double qf2_phys=1.0/9.0;
  djack_t Zv_e2_num=djvec_t(Spect_Spect+qf2_phys*e2*djvec_t(Spect_Spect_der_e2-deltam_cr*Spect_Spect_der_KA+deltam*Spect_Spect_der_MA))[T/2];
  djvec_t Zv_e2_den=Spect_Seq+qf2_phys*e2*djvec_t(Spect_Seq_der_e2-deltam_cr*Spect_Seq_der_KA+deltam*Spect_Seq_der_MA);
  djvec_t Zv_e2=Zv_e2_num/(-2.0*Zv_e2_den);
  Zv_e2.ave_err().write("plots/Zv_e2.xmg");
  djvec_t(Zv_e2-Zv).ave_err().write("plots/delta_Zv_e2.xmg");
  djvec_t((Zv_e2-Zv)/(Zv*qf2_phys)).ave_err().write("plots/Zv_e2_fact.xmg");
  
  double pred_Zv=-20.6178*alpha_em/(4*M_PI);
  double pred_Za=-15.7963*alpha_em/(4*M_PI);
  cout<<"Prediction from Davide, Zv: "<<pred_Zv<<endl;
  cout<<"Prediction from Davide, Za: "<<pred_Za<<endl;
  
  djack_t CZv_e2_num=djvec_t(Spect_Spect+qf2_phys*e2*djvec_t(Spect_Spect_der_e2-deltam_cr*Spect_Spect_der_KA+deltam*Spect_Spect_der_MA))[T/2];
  cout<<"CZv_e2_num: "<<CZv_e2_num.ave_err()<<endl;
  djvec_t CZv_e2_den=Spect_V_Seq+qf2_phys*e2*djvec_t(Spect_V_Seq_der_e2-deltam_cr*Spect_V_Seq_der_KA+deltam*Spect_V_Seq_der_MA);
  CZv_e2_den.ave_err().write("plots/CZv_e2_den.xmg");
  djvec_t CZv_e2=CZv_e2_num/(-2.0*CZv_e2_den);
  CZv_e2.ave_err().write("plots/CZv_e2.xmg");
  djvec_t((CZv_e2-CZv)/(CZv*qf2_phys)).ave_err().write("plots/CZv_e2_fact.xmg");
  
  return 0;
}
