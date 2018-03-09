#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/perens.hpp>

const int ODD=-1,EVN=1;

djvec_t perens_t::get_contraction(const string &combo,const string &ID,const size_t reim,const int tpar)
{
  const size_t T=L[0];
  
  const vector<string> cID={"V0P5","P5P5"};
  
  const size_t ncols=2;
  const string templ=dir_path+"/"+prop_hadr_path+"/%04zu/mes_contr_"+combo;
  const djvec_t data=read_conf_set_t(templ,conf_list,ncols,{0,1},T,SILENT);
  if(data.size()==0)
    {
      cout<<conf_list<<endl;
      CRASH("No file opened for template %s",templ.c_str());
    }
  
  const auto pos=find(cID.begin(),cID.end(),ID);
  if(pos==cID.end()) CRASH("Unknown %s",ID.c_str());
  
  //filter
  const size_t offset=ncols*distance(cID.begin(),pos)+reim;
  const size_t each=ncols*cID.size();
  const size_t base_nel=T;
  const size_t hw=data.size()/(base_nel*each);
  
  return vec_filter(data,gslice(base_nel*offset,{hw,T},{each*base_nel,1})).symmetrized(tpar);
}

djack_t perens_t::compute_deltam_cr(const size_t im)
{
  auto get=[im,this]
    (string tag_bw,string tag_fw,const string &ID,const size_t reim,const int tpar,const int rpar,const int rdiff)
    {
      djvec_t res(L[0]/2+1);
      res=0.0;
      for(size_t r=0;r<nr;r++)
	{
	  const string name="M"+to_string(im)+"_R"+to_string(r)+"_"+tag_bw+"_M"+to_string(im)+"_R"+to_string((r+rpar)%nr)+"_"+tag_fw;
	  const djvec_t contr=get_contraction(name,ID,reim,tpar)*((r==0)?1:rpar)/(1+abs(rpar));
	  contr.ave_err().write(dir_path+"/plots/"+ID+"_"+name+".xmg");
	  res+=contr;
	}
      
      return res;
    };
  
  const djvec_t P5P5_00=get("0","0","P5P5",RE,EVN,EVN,1);
  const djack_t m_P=constant_fit(effective_mass(P5P5_00),tmin,tmax,dir_path+"/plots/m_P_"+to_string(im)+".xmg");
  cout<<"M["<<im<<"]: "<<m_P<<endl;
  
  //measure mcrit according to eq.3 of hep-lat/0701012
  const djvec_t V0P5_00=get("0","0","V0P5",IM,ODD,ODD,1);
  const djvec_t m_cr_corr=forward_derivative(V0P5_00)/(2.0*P5P5_00);
  const djack_t m_cr=constant_fit(m_cr_corr,tmin,tmax,dir_path+"/plots/m_cr_"+to_string(im)+".xmg");
  
  if(pars::use_QED)
    {
      //load corrections
      const djvec_t V0P5_LL=get("F","F","V0P5",IM,ODD,ODD,1);
      const djvec_t V0P5_0M=get("0","FF","V0P5",IM,ODD,ODD,1);
      const djvec_t V0P5_M0=get("FF","0","V0P5",IM,ODD,ODD,1);
      const djvec_t V0P5_0T=get("0","T","V0P5",IM,ODD,ODD,1);
      const djvec_t V0P5_T0=get("T","0","V0P5",IM,ODD,ODD,1);
      //load the derivative wrt counterterm
      const djvec_t V0P5_0P=get("0","P","V0P5",RE,ODD,EVN,1);
      const djvec_t V0P5_P0=get("P","0","V0P5",RE,ODD,EVN,1);
      
      //build numerator
      const djvec_t num_deltam_cr_corr=V0P5_LL+V0P5_0M+V0P5_M0+V0P5_0T+V0P5_T0;
      //build denominator
      const djvec_t den_deltam_cr_corr=V0P5_P0-V0P5_0P;
      const djvec_t deltam_cr_corr=num_deltam_cr_corr/den_deltam_cr_corr;
      const djack_t deltam_cr=constant_fit(deltam_cr_corr,tmin,tmax,dir_path+"/plots/deltam_cr.xmg");
      
      return deltam_cr;
    }
  else return djack_t{};
}

djack_t perens_t::compute_meson_mass(const size_t im1,const size_t im2)
{
  djvec_t P5P5_corr(L[0]/2+1);
  P5P5_corr=0.0;
  for(size_t r=0;r<nr;r++)
    {
      string name="M"+to_string(im1)+"_R"+to_string(r)+"_0_M"+to_string(im2)+"_R"+to_string(r)+"_0";
      djvec_t contr=get_contraction(name,"P5P5",RE,EVN);
      P5P5_corr+=contr;
    }
  P5P5_corr/=nr;
  
  const djack_t m_P=constant_fit(effective_mass(P5P5_corr),tmin,tmax,dir_path+"/plots/m_P_"+to_string(im1)+"_"+to_string(im1)+".xmg");
  cout<<"M["<<im1<<","<<im2<<"]: "<<m_P<<endl;
  
  return m_P;
}
