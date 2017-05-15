#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=48;
const size_t nconfs=10;
const range_t conf_set={700,10,700+35*(nconfs-1)};

//! return the whole correlator
djvec_t get(const string &combo,const string &ID,const size_t reim,const size_t sym)
{
  const vector<string> cID=
    {"S0P5","V1P5","V2P5","V3P5","V0P5","P5P5","A1P5","A2P5","A3P5","A0P5","T1P5","T2P5","T3P5","B1P5","B2P5","B3P5","P5S0","P5V1","P5V2","P5V3","P5V0",
     "P5P5","P5A1","P5A2","P5A3","P5A0","P5T1","P5T2","P5T3","P5B1","P5B2","P5B3","P5P5"};
 
 const size_t ncols=2;
 const string templ="%04d/mes_contr_"+combo;
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
 res.ave_err().write("plots/"+combo+".xmg");
 
 //self-symmetrize
 res.symmetrize(sym);
 
 //write symmetrized
 res.ave_err().write("plots/"+combo+"_symm.xmg");
 
 return res;
}

int main()
{
  const double m=0.02363,mp=0.02373,mm=0.02353;
  const double k=0.163255,kp=0.163355,km=0.163155;
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  
  set_njacks(nconfs);
  
  djvec_t t1=get("Source_Seq","S0P5",RE,1);
  cout<<t1[0].ave_err()<<endl;
  
  djvec_t t2=get("Spect_Spect","P5P5",RE,1);
  cout<<t2[T/2].ave_err()<<endl;
  
  get("Spect_Spect","P5P5",RE,1);
  get("SpectP_Spect","P5P5",RE,1);
  
  //2pts
  djvec_t Spect_Spect=get("Spect_Spect","P5P5",RE,1);
  djack_t M,Z2;
  two_pts_fit(Z2,M,Spect_Spect,T/2,12,T/2,"plots/Spect_Spect_effmass.xmg");
  cout<<"Mass: "<<M.ave_err()<<endl;
  
  //3pts
  djvec_t Spect_Seq=get("Spect_Seq","V0P5",RE,-1);
  djvec_t Spect_CV_Seq=get("Spect_V_Seq","S0P5",RE,-1);
  
  //3pts derivative w.r.t to e2
  djvec_t SpectP_SeqP=get("SpectP_SeqP","V0P5",RE,-1);
  djvec_t SpectM_SeqM=get("SpectM_SeqM","V0P5",RE,-1);
  djvec_t Spect_Seq_der_P=djvec_t(SpectP_SeqP-Spect_Seq)/(sqr(0.25)*e2);
  Spect_Seq_der_P.ave_err().write("plots/Spect_Seq_der_P.xmg");
  djvec_t Spect_Seq_der_M=djvec_t(SpectM_SeqM-Spect_Seq)/(sqr(0.25)*e2);
  Spect_Seq_der_M.ave_err().write("plots/Spect_Seq_der_M.xmg");
  djvec_t Spect_Seq_der_=djvec_t(Spect_Seq_der_P+Spect_Seq_der_M)/2;
  Spect_Seq_der_.ave_err().write("plots/Spect_Seq_der_.xmg");
  
  //3pts derivative w.r.t mass
  djvec_t SpectMAP_SeqMAP=get("SpectMAP_SeqMAP","V0P5",RE,-1);
  djvec_t SpectMAM_SeqMAM=get("SpectMAM_SeqMAM","V0P5",RE,-1);
  djvec_t Spect_Seq_der_MAP=djvec_t(SpectMAP_SeqMAP-Spect_Seq)/(mp-m);
  Spect_Seq_der_MAP.ave_err().write("plots/Spect_Seq_der_MAP.xmg");
  djvec_t Spect_Seq_der_MAM=djvec_t(SpectMAM_SeqMAM-Spect_Seq)/(mm-m);
  Spect_Seq_der_MAM.ave_err().write("plots/Spect_Seq_der_MAM.xmg");
  djvec_t Spect_Seq_der_MA=djvec_t(Spect_Seq_der_MAP+Spect_Seq_der_MAM)/2;
  Spect_Seq_der_MA.ave_err().write("plots/Spect_Seq_der_MA.xmg");
  
  //3pts derivative w.r.t critical mass
  djvec_t SpectKAP_SeqKAP=get("SpectKAP_SeqKAP","V0P5",RE,-1);
  djvec_t SpectKAM_SeqKAM=get("SpectKAM_SeqKAM","V0P5",RE,-1);
  djvec_t Spect_Seq_der_KAP=djvec_t(SpectKAP_SeqKAP-Spect_Seq)/(mcp-mc);
  Spect_Seq_der_KAP.ave_err().write("plots/Spect_Seq_der_KAP.xmg");
  djvec_t Spect_Seq_der_KAM=djvec_t(SpectKAM_SeqKAM-Spect_Seq)/(mcm-mc);
  Spect_Seq_der_KAM.ave_err().write("plots/Spect_Seq_der_KAM.xmg");
  djvec_t Spect_Seq_der_KA=djvec_t(Spect_Seq_der_KAP+Spect_Seq_der_KAM)/2;
  Spect_Seq_der_KA.ave_err().write("plots/Spect_Seq_der_KA.xmg");
  
  //Zv
  djvec_t(Spect_Spect[T/2]/(-2.0*Spect_Seq)).ave_err().write("plots/Spect_Seq_Zv.xmg");
  djvec_t(Spect_Spect[T/2]/(-2.0*Spect_CV_Seq)).ave_err().write("plots/Spect_Seq_CZv.xmg");
  
  //3pts derivative w.r.t to e2
  djvec_t SpectP_CV_SeqP=get("SpectP_V_SeqP","V0P5",RE,1);
  djvec_t SpectM_CV_SeqM=get("SpectM_V_SeqM","V0P5",RE,1);
  djvec_t Spect_CV_Seq_der_P=djvec_t(SpectP_CV_SeqP-Spect_CV_Seq)/(sqr(0.25)*e2);
  Spect_CV_Seq_der_P.ave_err().write("plots/Spect_CV_Seq_der_P.xmg");
  djvec_t Spect_CV_Seq_der_M=djvec_t(SpectM_CV_SeqM-Spect_CV_Seq)/(sqr(0.25)*e2);
  Spect_CV_Seq_der_M.ave_err().write("plots/Spect_CV_Seq_der_M.xmg");
  djvec_t Spect_CV_Seq_der_=djvec_t(Spect_CV_Seq_der_P+Spect_CV_Seq_der_M)/2;
  Spect_CV_Seq_der_.ave_err().write("plots/Spect_CV_Seq_der_.xmg");
  
  return 0;
}
