#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/perens.hpp>

djvec_t perens_t::get_contraction_by_name(const string &suffix,const string &bil_name,const dcompl_t &c_coeff,const int tpar)
{
  const size_t T=L[0];
  
  const vector<string> cID={"V0P5","P5P5"};
  
  const size_t ncols=2;
  const string templ=dir_path+"/"+prop_hadr_path+"/%04zu/mes_contr_"+suffix;
  const djvec_t data=read_conf_set_t(templ,conf_list,ncols,{0,1},T,SILENT);
  if(data.size()==0)
    {
      cout<<conf_list<<endl;
      CRASH("No file opened for template %s",templ.c_str());
    }
  
  const auto pos=find(cID.begin(),cID.end(),bil_name);
  if(pos==cID.end()) CRASH("Unknown %s",bil_name.c_str());
  
  if(real(c_coeff)!=0 and imag(c_coeff)!=0) CRASH("Don't know what to do");
  
  //Find whether to take real or imaginary, and the coefficient
  size_t reim;
  double coeff;
  if(real(c_coeff)!=0)
    {
      reim=0;
      coeff=real(c_coeff);
    }
  else
    {
      reim=1;
      coeff=-imag(c_coeff);
    }
  
  //filter
  const size_t offset=ncols*distance(cID.begin(),pos)+reim;
  const size_t each=ncols*cID.size();
  const size_t base_nel=T;
  const size_t hw=data.size()/(base_nel*each);
  
  const djvec_t out=coeff*vec_filter(data,gslice(base_nel*offset,{hw,T},{each*base_nel,1}));
  
  if(tpar) return out.symmetrized(tpar);
  else     return out;
}

djvec_t perens_t::get_contraction(const int imbw,qprop::ins kbw,const int imfw,qprop::ins kfw,const string &ID,const size_t ext_reim,const int tpar,const size_t rfw,const int rdiff)
{
  const string tag_bw=qprop::ins_tag[kbw];
  const string tag_fw=qprop::ins_tag[kfw];
  
  const size_t rbw=(rfw+rdiff)%nr;
  
  //Compute the coefficient
  dcompl_t c_coeff=
    conj(coeff_to_read(kbw,rbw))*
    coeff_to_read(kfw,rfw);
  
  const string m_tag_bw=(imbw==(int)nm)?"sea":to_string(imbw);
  const string m_tag_fw=(imfw==(int)nm)?"sea":to_string(imfw);
  
  //Include -i if asking the imaginary part
  if(ext_reim==1) c_coeff*=dcompl_t(0.0,-1.0);
  
  const string name="M"+m_tag_bw+"_R"+to_string(rbw)+"_"+tag_bw+"_M"+m_tag_fw+"_R"+to_string(rfw)+"_"+tag_fw;
  const djvec_t res=get_contraction_by_name(name,ID,c_coeff,tpar);
  res.ave_err().write(dir_path+"/plots/"+ID+"_"+name+".xmg");
  
  return res;
}

djack_t perens_t::compute_meson_mass(const string& m1_tag,const string& m2_tag)
{
  djvec_t P5P5_corr(L[0]/2+1);
  P5P5_corr=0.0;
  for(size_t r=0;r<nr;r++)
    {
      string name="M"+m1_tag+"_R"+to_string(r)+"_0_M"+m2_tag+"_R"+to_string(r)+"_0";
      djvec_t contr=get_contraction_by_name(name,"P5P5",1.0,EVN);
      P5P5_corr+=contr;
    }
  P5P5_corr/=nr;
  
  const djack_t m_P=constant_fit(effective_mass(P5P5_corr),tmin,tmax,dir_path+"/plots/m_P_"+m1_tag+"_"+m1_tag+".xmg");
  cout<<"M["<<m1_tag<<","<<m2_tag<<"]: "<<m_P<<endl;
  
  return m_P;
}

djack_t perens_t::compute_meson_mass_QED(const size_t im1,const size_t im2)
{
  djvec_t P5P5_corr(L[0]/2+1);
  djvec_t P5P5_QED_corr(L[0]/2+1);
  P5P5_corr=P5P5_QED_corr=0.0;
  for(size_t r=0;r<nr;r++)
    {
      using namespace qprop;
      const size_t rdiff=0;
      auto cg=[r,this](size_t im1,size_t im2,ins tfw,ins tbw){return get_contraction(im1,tfw,im2,tbw,"P5P5",RE,EVN,r,rdiff);};
      P5P5_corr+=cg(im1,im2,LO,LO);
      switch(pars::use_QED)
	{
	case 0:
	  break;
	case 1:
	  CRASH("Not implemented");
	  break;
	case 2:
	  P5P5_QED_corr+=cg(im1,im2,F,F)+cg(im1,im2,LO,QED)+cg(im2,im1,LO,QED); //hack
	  break;
	}
    }
  P5P5_corr/=nr;
  P5P5_QED_corr/=nr;
  
  const djack_t m_P=constant_fit(effective_mass(P5P5_corr),tmin,tmax,dir_path+"/plots/m_P_"+to_string(im1)+"_"+to_string(im2)+".xmg");
  const djack_t dm_P=-constant_fit(effective_slope(djvec_t(P5P5_QED_corr/P5P5_corr),effective_mass(P5P5_corr),L[0]/2),tmin,tmax,dir_path+"/plots/dm_P_"+to_string(im1)+"_"+to_string(im2)+".xmg");
  
  cout<<"dM["<<im1<<","<<im2<<"]: "<<dm_P<<endl;
  
  return dm_P;
}

djack_t perens_t::compute_mPCAC(const size_t im)
{
  djvec_t P5P5_LO(L[0]/2+1);
  djvec_t V0P5_LO(L[0]/2+1);
  djvec_t P5P5_QED(L[0]/2+1);
  djvec_t V0P5_QED(L[0]/2+1);
  P5P5_LO=V0P5_LO=P5P5_QED=V0P5_QED=0.0;
  const size_t rdiff=0;
  for(size_t r=0;r<nr;r++)
    {
      using namespace qprop;
      P5P5_LO+=get_contraction(im,LO,im,LO,"P5P5",RE,EVN,r,rdiff);
      V0P5_LO+=tau3[r]*get_contraction(im,LO,im,LO,"V0P5",IM,ODD,r,rdiff);
      if(pars::use_QED==2 and im!=nm)
	{
	  P5P5_QED+=2*get_contraction(im,LO,im,QED,"P5P5",RE,ODD,r,rdiff);
	  V0P5_QED+=2*tau3[r]*get_contraction(im,LO,im,QED,"V0P5",IM,EVN,r,rdiff);
	}
    }
  P5P5_LO/=nr;
  V0P5_LO/=nr;
  P5P5_QED/=nr;
  V0P5_QED/=nr;
  
  const string m_tag=to_string(im);
  const djvec_t mPCAC_LO=forward_derivative(V0P5_LO)/(2.0*P5P5_LO);
  const djvec_t mPCAC_LO_symm=symmetric_derivative(V0P5_LO)/(2.0*P5P5_LO);
  const djack_t mPCAC=constant_fit(mPCAC_LO,tmin,tmax,dir_path+"/plots/mPCAC_"+m_tag+".xmg");
  const djack_t mPCAC_symm=constant_fit(mPCAC_LO_symm,tmin,tmax,dir_path+"/plots/mPCAC_symm"+m_tag+".xmg");
  cout<<"MPCAC["<<im<<"]:      "<<smart_print(mPCAC)<<endl;
  cout<<"MPCAC["<<im<<"] symm: "<<smart_print(mPCAC_symm)<<endl;
  
  if(pars::use_QED==2 and im!=nm)
    {
      const djvec_t a=djvec_t(symmetric_derivative(V0P5_QED)/P5P5_LO-symmetric_derivative(V0P5_LO)/sqr(P5P5_LO)*P5P5_QED);
      const djack_t dM_pcac=constant_fit(a,tmin,tmax,dir_path+"/plots/var_m_pcac_m"+m_tag+"_P5P5.xmg");
      cout<<"dM_pcac, M["<<im<<"]: "<<smart_print(dM_pcac)<<endl;
    }
  
  return mPCAC_symm;
}
