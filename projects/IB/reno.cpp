#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=48;
const size_t nconfs=10;
const range_t conf_set={700,10,700+35*(nconfs-1)};

//! return the whole correlator
djvec_t get(const string &combo,const string &ID,const size_t reim)
{
const vector<string> cID=
  {"S0P5","V1P5","V2P5","V3P5","V0P5","P5P5","A1P5","A2P5","A3P5","A0P5","T1P5","T2P5","T3P5","B1P5","B2P5","B3P5","P5S0","P5V1","P5V2","P5V3","P5V0",
   "P5P5","P5A1","P5A2","P5A3","P5A0","P5T1","P5T2","P5T3","P5B1","P5B2","P5B3","P5P5"};
 
 const size_t ncols=2;
 const djvec_t data=read_conf_set_t("%04d/mes_contr_"+combo,conf_set,ncols,{0,1},T,SILENT);
 if(data.size()==0) CRASH("No file opened");
 
 const auto pos=find(cID.begin(),cID.end(),ID);
 if(pos==cID.end()) CRASH("Unknown %s",ID.c_str());
 
 const size_t offset=ncols*distance(cID.begin(),pos)+reim;
 const size_t each=ncols*cID.size();
 const size_t base_nel=T;
 const size_t hw=data.size()/(base_nel*each);
 
 return vec_filter(data,gslice(base_nel*offset,{hw,T},{each*base_nel,1}));
}

int main()
{
  const double m=0.02363,mp=0.02373,mm=0.02353;
  const double k=0.163255,kp=0.163355,km=0.163155;
  const double mc=1/(2*k),mcp=1/(2*kp),mcm=1/(2*km);
  
  set_njacks(nconfs);
  
  djvec_t t1=get("Source_Seq","S0P5",RE);
  cout<<t1[0].ave_err()<<endl;
  
  djvec_t t2=get("Spect_Spect","P5P5",RE);
  cout<<t2[T/2].ave_err()<<endl;
  
  get("Spect_Spect","P5P5",RE).ave_err().write("Spect_Spect.xmg");
  get("SpectP_Spect","P5P5",RE).ave_err().write("SpectP_Spect.xmg");
  
  djvec_t Spect_Spect=get("Spect_Spect","P5P5",RE);
  Spect_Spect.ave_err().write("Spect_Spect.xmg");
  
  //derivative w.r.t mass
  djvec_t SpectP_Spect=get("SpectP_Spect","P5P5",RE);
  SpectP_Spect.ave_err().write("SpectP_Spect.xmg");
  djvec_t SpectM_Spect=get("SpectM_Spect","P5P5",RE);
  SpectM_Spect.ave_err().write("SpectM_Spect.xmg");
  djvec_t Spect_Spect_der_P=djvec_t(SpectP_Spect-Spect_Spect);
  Spect_Spect_der_P.ave_err().write("Spect_Spect_der_P.xmg");
  djvec_t Spect_Spect_der_M=djvec_t(SpectM_Spect-Spect_Spect);
  Spect_Spect_der_M.ave_err().write("Spect_Spect_der_M.xmg");
  djvec_t Spect_Spect_der_=djvec_t(Spect_Spect_der_P+Spect_Spect_der_M)/2;
  Spect_Spect_der_.ave_err().write("Spect_Spect_der_.xmg");
  
  //derivative w.r.t mass
  djvec_t SpectMAP_Spect=get("SpectMAP_Spect","P5P5",RE);
  SpectMAP_Spect.ave_err().write("SpectMAP_Spect.xmg");
  djvec_t SpectMAM_Spect=get("SpectMAM_Spect","P5P5",RE);
  SpectMAM_Spect.ave_err().write("SpectMAM_Spect.xmg");
  djvec_t Spect_Spect_der_MAP=djvec_t(SpectMAP_Spect-Spect_Spect)/(mp-m);
  Spect_Spect_der_MAP.ave_err().write("Spect_Spect_der_MAP.xmg");
  djvec_t Spect_Spect_der_MAM=djvec_t(SpectMAM_Spect-Spect_Spect)/(mm-m);
  Spect_Spect_der_MAM.ave_err().write("Spect_Spect_der_MAM.xmg");
  djvec_t Spect_Spect_der_MA=djvec_t(Spect_Spect_der_MAP+Spect_Spect_der_MAM)/2;
  Spect_Spect_der_MA.ave_err().write("Spect_Spect_der_MA.xmg");
  
  //derivative w.r.t  critical mass
  djvec_t SpectKAP_Spect=get("SpectKAP_Spect","P5P5",RE);
  SpectKAP_Spect.ave_err().write("SpectKAP_Spect.xmg");
  djvec_t SpectKAM_Spect=get("SpectKAM_Spect","P5P5",RE);
  SpectKAM_Spect.ave_err().write("SpectKAM_Spect.xmg");
  djvec_t Spect_Spect_der_KAP=djvec_t(SpectKAP_Spect-Spect_Spect)/(mcp-mc);
  Spect_Spect_der_KAP.ave_err().write("Spect_Spect_der_KAP.xmg");
  djvec_t Spect_Spect_der_KAM=djvec_t(SpectKAM_Spect-Spect_Spect)/(mcm-mc);
  Spect_Spect_der_KAM.ave_err().write("Spect_Spect_der_KAM.xmg");
  djvec_t Spect_Spect_der_KA=djvec_t(Spect_Spect_der_KAP+Spect_Spect_der_KAM)/2;
  Spect_Spect_der_KA.ave_err().write("Spect_Spect_der_KA.xmg");
  
  djvec_t Spect_Seq_P5V0P5=get("Spect_Seq","V0P5",RE).symmetrized(-1);
  Spect_Seq_P5V0P5.ave_err().write("Spect_Seq_P5V0P5.xmg");
  djvec_t(Spect_Spect[T/2]/(-2.0*Spect_Seq_P5V0P5)).ave_err().write("Spect_Seq_Zv.xmg");
  
  djvec_t Spect_Seq_P5VCP5=get("Spect_Seq","S0P5",RE).symmetrized(1);
  Spect_Seq_P5VCP5.ave_err().write("Spect_Seq_P5VCP5.xmg");
  djvec_t(2*Spect_Spect[T/2]/(Spect_Seq_P5VCP5)).ave_err().write("Spect_Seq_CZv.xmg");
  
  return 0;
}
