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

void set_mr_ind(size_t nm,size_t nr)
{
  mr_ind.set_ranges({{"im",nm},{"r",nr}});
  nmr=mr_ind.max();
}

string get_prop_tag(size_t im,size_t ir,const string &ins)
{return combine("S_M%zu_R%zu_%s",im,ir,ins.c_str());}

void set_jprops(bool set_QED)
{
  cout<<"Setting all "<<nmr<<" jprops"<<endl;
  
  // jprop_0.resize(nmr);
  // if(set_QED)
  //   for(auto &o : {&jprop_2,&jprop_P,&jprop_S})
  //     o->resize(nmr);
}

void clusterize_all_mr_INS_props(vjprop_t &jprop,size_t clust_size)
{
#pragma omp parallel for
  for(size_t imr=0;imr<nmr;imr++)
    clusterize(jprop[imr],clust_size);
}

void clusterize_all_mr_props(bool use_QED,size_t clust_size)
{
  cout<<"Clusterizing all props, clust_size="<<clust_size<<endl;
  
  // clusterize_all_mr_INS_props(jprop_0,clust_size);
  // if(use_QED)
  //   for(auto &p : {&jprop_2,&jprop_P,&jprop_S})
  //     clusterize_all_mr_INS_props(*p,clust_size);
}

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

