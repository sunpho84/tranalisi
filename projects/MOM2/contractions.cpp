#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/contractions.hpp>
#include <MOM2/perens.hpp>

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

djack_t perens_t::compute_meson_mass(const string& m1_tag,const string& m2_tag)
{
  djvec_t P5P5_corr(L[0]/2+1);
  P5P5_corr=0.0;
  for(size_t r=0;r<nr;r++)
    {
      string name="M"+m1_tag+"_R"+to_string(r)+"_0_M"+m2_tag+"_R"+to_string(r)+"_0";
      djvec_t contr=get_contraction(name,"P5P5",RE,EVN);
      P5P5_corr+=contr;
    }
  P5P5_corr/=nr;
  
  const djack_t m_P=constant_fit(effective_mass(P5P5_corr),tmin,tmax,dir_path+"/plots/m_P_"+m1_tag+"_"+m1_tag+".xmg");
  cout<<"M["<<m1_tag<<","<<m2_tag<<"]: "<<m_P<<endl;
  
  return m_P;
}