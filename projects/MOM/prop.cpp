#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PROP
 #include <prop.hpp>

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

void read_prop(prop_t &prop,raw_file_t &file,const dcompl_t &fact)
{
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    prop(isc(is_si,ic_si),isc(is_so,ic_so))=c*fact;
	  }
}

string get_prop_tag(size_t im,size_t ir,const string &ins)
{return combine("S_M%zu_R%zu_%s",im,ir,ins.c_str());}

vjprop_t get_all_mr_props_inv(const vjprop_t &jprop)
{
  cout<<"Inverting all props"<<endl;
  
  vjprop_t jprop_inv(jprop.size());
  
#pragma omp parallel for
  for(size_t imr=0;imr<nmr;imr++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      put_into_jackknife(jprop_inv[imr],get_from_jackknife(jprop[imr],ijack).inverse(),ijack);
  
  return jprop_inv;
}

void build_all_mr_jackknifed_props(vector<jm_r_mom_props_t> &jprops,const vector<m_r_mom_conf_props_t> &props,bool set_QED,const index_t &im_r_ijack_ind,const index_t &im_r_ind)
{
  for(size_t i_im_r_ijack=0;i_im_r_ijack<im_r_ijack_ind.max();i_im_r_ijack++)
    {
      vector<size_t> im_r_ijack=im_r_ijack_ind(i_im_r_ijack);
      size_t im=im_r_ijack[0],r=im_r_ijack[1],ijack=im_r_ijack[2];
      size_t im_r=im_r_ind({im,r});
      
      jm_r_mom_props_t &j=jprops[im_r];
      const m_r_mom_conf_props_t &p=props[i_im_r_ijack];
      
      add_to_cluster(j.jprop_0,p.prop_0,ijack);
      if(set_QED)
	for(auto &jp_p : vector<pair<jprop_t*,const prop_t*>>({{&j.jprop_2,&p.prop_FF},{&j.jprop_2,&p.prop_T},{&j.jprop_P,&p.prop_P},{&j.jprop_S,&p.prop_S}}))
	  add_to_cluster(*jp_p.first,*jp_p.second,ijack);
    }
}
