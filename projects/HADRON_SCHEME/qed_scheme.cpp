#include <tranalisi.hpp>

size_t T;
range_t file_range;

djvec_t read(const string& name,const size_t which,const bool reIm,const double sign,const int parity)
{
  const djvec_t data=read_conf_set_t("out/%04d/"+name,file_range,2,{reIm},T,true);
  const size_t base_nel=T;
  const size_t offset=which;
  const size_t each=2;
  const size_t hw=data.size()/(base_nel*each);
  
  return sign*vec_filter(data,gslice(base_nel*offset,{hw,T},{each*T,1})).symmetrized(parity);
}

djvec_t load_bar(const string& tag,const size_t iq,const int iproj,const int iWick,const bool reIm,const double sign,const int parity)
{
  string name=combine("bar_alt_contr_SM_q%zu_q%zu_q%zu_%s_proj_%d_Wick_%d",iq,iq,iq,tag.c_str(),iproj,iWick);
  
  return -read(name,0,reIm,sign,parity);
}

djvec_t load_mes(const string& plotName,const size_t which,const char *tag_bw,const size_t iq_fw,const char *tag_fw,const bool reIm,const double sign,const int parity)
{
  string name=combine("mes_contr_M0_R0_%s_M%zu_R0_%s",tag_bw,iq_fw,tag_fw);
  
  const djvec_t out=read(name,which,reIm,sign,parity);
  
  out.ave_err().write("plots/"+plotName+".xmg");
  
  return out;
}

int main()
{
  set_njacks(15);
  file_range.start=700;
  file_range.each=10;
  file_range.end=2190;
  T=48;
  
  const int EVN=+1,ODD=-1,V0P5=0,P5P5=1;
  const int qfw=1,ifw=0;
  const int qbw=0;//ibw=0 always
  
  const djvec_t P5P5_00=load_mes("P5P5_00",P5P5,"0",ifw,"0",RE,+1,EVN);
  
  const djvec_t P5P5_QQ=load_mes("P5P5_QQ",P5P5,"F",ifw,"F",RE,+1,EVN);
  const djvec_t P5P5_0Q=load_mes("P5P5_0Q",P5P5,"0",ifw,"QED",RE,+1,EVN);
  const djvec_t P5P5_Q0=load_mes("P5P5_Q0",P5P5,"QED",ifw,"0",RE,+1,EVN);
  const djvec_t P5P5_Q=qbw*qbw*P5P5_Q0+qfw*qfw*P5P5_0Q+qfw*qbw*P5P5_QQ;
  
  const djvec_t P5P5_0P=load_mes("P5P5_0P",P5P5,"0",ifw,"PSE",IM,-1,EVN);
  const djvec_t P5P5_P0=load_mes("P5P5_P0",P5P5,"PSE",ifw,"0",IM,+1,EVN);
  const djvec_t P5P5_P=qbw*P5P5_P0+qfw*P5P5_0P;
  
  const djvec_t P5P5_0S=load_mes("P5P5_0S",P5P5,"0",ifw,"MASS",RE,+1,EVN);
  const djvec_t P5P5_S0=load_mes("P5P5_S0",P5P5,"MASS",ifw,"0",RE,-1,EVN);
  
  const djvec_t V0P5_00=load_mes("V0P5_00",V0P5,"0",ifw,"0",IM,+1,ODD);
  
  const djvec_t V0P5_QQ=load_mes("V0P5_QQ",V0P5,"F",ifw,"F",IM,+1,ODD);
  const djvec_t V0P5_0Q=load_mes("V0P5_0Q",V0P5,"0",ifw,"QED",IM,+1,ODD);
  const djvec_t V0P5_Q0=load_mes("V0P5_Q0",V0P5,"QED",ifw,"0",IM,+1,ODD);
  const djvec_t V0P5_Q=qbw*qbw*V0P5_Q0+qfw*qfw*V0P5_0Q+qfw*qbw*V0P5_QQ;
  
  const djvec_t V0P5_0S=load_mes("V0P5_0S",V0P5,"0",ifw,"MASS",IM,+1,ODD);
  const djvec_t V0P5_S0=load_mes("V0P5_S0",V0P5,"MASS",ifw,"0",IM,+1,ODD);

  const djvec_t V0P5_0P=load_mes("V0P5_0P",V0P5,"0",ifw,"PSE",RE,+1,ODD);
  const djvec_t V0P5_P0=load_mes("V0P5_P0",V0P5,"PSE",ifw,"0",RE,-1,ODD);
  const djvec_t V0P5_P=qbw*V0P5_P0+qfw*V0P5_0P;
  
  const djvec_t mcr_corr=forward_derivative(V0P5_00)/P5P5_00;
  const djvec_t dmcr_corr_P=-forward_derivative(V0P5_P)/P5P5_00+P5P5_P*forward_derivative(V0P5_00)/sqr(P5P5_00);
  const djvec_t dmcr_corr_Q=-forward_derivative(V0P5_Q)/P5P5_00+P5P5_Q*forward_derivative(V0P5_00)/sqr(P5P5_00);
  const djvec_t dP_corr=forward_derivative(V0P5_Q)/forward_derivative(V0P5_P);
  const djvec_t dP_full_corr=-dmcr_corr_Q/dmcr_corr_P;
  
  const djack_t mcr=constant_fit(mcr_corr,10,15,"plots/mcr_corr.xmg");
  const djack_t dP=constant_fit(dP_corr,10,15,"plots/dP_corr.xmg");
  const djack_t dP_full=constant_fit(dP_full_corr,10,15,"plots/dP_full_corr.xmg");
  
  cout<<"dP: "<<smart_print(dP)<<endl;
  cout<<"dP_full: "<<smart_print(dP_full)<<endl;
  
  return 0;
}

