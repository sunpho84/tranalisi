#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <fit.hpp>
#include <oper.hpp>

#include <contractions.hpp>
#include <geometry.hpp>

const int ODD=-1,EVN=1;

djvec_t get_contraction(const string &combo,const string &ID,vector<size_t> &conf_list,const size_t reim,const int tpar)
{
  const size_t T=L[0];
  
  const vector<string> cID={"V0P5","P5P5"};
  
  const size_t ncols=2;
  const string templ="out/%04zu/mes_contr_"+combo;
  const djvec_t data=read_conf_set_t(templ,conf_list,ncols,{0,1},T,SILENT);
  if(data.size()==0) CRASH("No file opened for template %s",templ.c_str());
  
  const auto pos=find(cID.begin(),cID.end(),ID);
  if(pos==cID.end()) CRASH("Unknown %s",ID.c_str());
  
  //filter
  const size_t offset=ncols*distance(cID.begin(),pos)+reim;
  const size_t each=ncols*cID.size();
  const size_t base_nel=T;
  const size_t hw=data.size()/(base_nel*each);
  
  return vec_filter(data,gslice(base_nel*offset,{hw,T},{each*base_nel,1})).symmetrized(tpar);
}

djack_t compute_deltam_cr(vector<size_t> &conf_list,const size_t tmin,const size_t tmax,const size_t im)
{
  auto get=[&conf_list,im]
    (string tag_bw,string tag_fw,const string &ID,const size_t reim,const int tpar)
    {
      string name="M"+to_string(im)+"_R0_"+tag_bw+"_M"+to_string(im)+"_R0_"+tag_fw;
      djvec_t res=get_contraction(name,ID,conf_list,reim,tpar);
      
      res.ave_err().write("plots/"+name+".xmg");
      
      return res;
    };
  
  djvec_t P5P5_00=get("0","0","P5P5",RE,EVN);
  
  //measure mcrit according to eq.3 of hep-lat/0701012
  djvec_t V0P5_00=get("0","0","V0P5",IM,ODD);
  djvec_t m_cr_corr=forward_derivative(V0P5_00)/(2.0*P5P5_00);
  djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,"plots/m_cr.xmg");
  
  //load corrections
  djvec_t V0P5_LL=get("F","F","V0P5",IM,ODD);
  djvec_t V0P5_0M=get("0","FF","V0P5",IM,ODD);
  djvec_t V0P5_M0=get("FF","0","V0P5",IM,ODD);
  djvec_t V0P5_0T=get("0","T","V0P5",IM,ODD);
  djvec_t V0P5_T0=get("T","0","V0P5",IM,ODD);
  //load the derivative wrt counterterm
  djvec_t V0P5_0P=get("0","P","V0P5",RE,ODD);
  djvec_t V0P5_P0=get("P","0","V0P5",RE,ODD);
  
  //build numerator
  djvec_t num_deltam_cr_corr=V0P5_LL+V0P5_0M+V0P5_M0+V0P5_0T+V0P5_T0+V0P5_0P+V0P5_P0;
  //build denominator
  djvec_t den_deltam_cr_corr=V0P5_P0-V0P5_0P;
  djvec_t deltam_cr_corr=num_deltam_cr_corr/den_deltam_cr_corr;
  djack_t deltam_cr=constant_fit(deltam_cr_corr,tmin,tmax,"plots/deltam_cr.xmg");
  
  return deltam_cr;
}
