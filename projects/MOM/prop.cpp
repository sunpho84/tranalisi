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

prop_t read_prop(raw_file_t &file)
{
  prop_t prop;
  
  //! source file
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t is_si=0;is_si<NSPIN;is_si++)
	for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	  {
	    dcompl_t c;
	    file.bin_read(c);
	    prop(isc(is_si,ic_si),isc(is_so,ic_so))=c;
	  }
  
  return prop;
}

void set_mr_ind(size_t nm,size_t nr)
{
  mr_ind.set_ranges({{"im",nm},{"r",nr}});
  nmr=mr_ind.max();
}

string get_prop_tag(size_t im,size_t ir,const string &ins)
{return combine("S_M%zu_R%zu_%s",im,ir,ins.c_str());}

void set_conf_props(bool set_QED)
{
  cout<<"Setting all "<<nmr<<" conf props"<<endl;
  
  //resize all props
  mom_prop_0.resize(nmr);
  if(set_QED)
    for(auto &mom_prop : {&mom_prop_FF,&mom_prop_F,&mom_prop_T,&mom_prop_P,&mom_prop_S})
      mom_prop->resize(nmr);
}

void read_all_mr_INS_props(vprop_t &mom_prop,map<string,vector<raw_file_t>> &map_files,const string &ins,size_t iconf_hit)
{
  if(map_files.find(ins)==map_files.end()) CRASH("Unable to find prop with insertion kind: %s",ins.c_str());
  for(size_t imr=0;imr<nmr;imr++) mom_prop[imr]=read_prop(map_files[ins][imr]);
}

void read_all_mr_props(bool read_QED,map<string,vector<raw_file_t>> &map_files,size_t iconf_hit)
{
  read_all_mr_INS_props(mom_prop_0,map_files,"0",iconf_hit);
  
  if(read_QED)
    {
      for(auto &o : vector<pair<vprop_t*,string>>({{&mom_prop_FF,"FF"},{&mom_prop_F,"F"},{&mom_prop_T,"T"},{&mom_prop_P,"P"},{&mom_prop_S,"S"}}))
	read_all_mr_INS_props(*o.first,map_files,o.second,iconf_hit);
      
      dcompl_t fact_P(0.0,-1.0);
      dcompl_t fact_S(-1.0,0.0);
      
      //put factors
      for(auto &pmr_f : vector<pair<vprop_t*,dcompl_t>>({{&mom_prop_P,fact_P},{&mom_prop_S,fact_S}}))
	for(auto &p : *pmr_f.first) //all mr
	  p*=pmr_f.second;
    }
}

void set_jprops(bool set_QED)
{
  cout<<"Setting all "<<nmr<<" jprops"<<endl;
  
  jprop_0.resize(nmr);
  if(set_QED)
    for(auto &o : {&jprop_2,&jprop_P,&jprop_S})
      o->resize(nmr);
}

void build_all_mr_jackknifed_INS_props(vjprop_t &out,const vprop_t &in,size_t ijack)
{
  for(size_t imr=0;imr<nmr;imr++)
    add_to_cluster(out[imr],in[imr],ijack);
}

void build_all_mr_jackknifed_props(bool set_QED,size_t ijack)
{
  build_all_mr_jackknifed_INS_props(jprop_0,mom_prop_0,ijack);
  if(set_QED)
    for(auto &jp_p : vector<pair<vjprop_t*,vprop_t*>>({{&jprop_2,&mom_prop_FF},{&jprop_2,&mom_prop_T},{&jprop_P,&mom_prop_P},{&jprop_S,&mom_prop_S}}))
  build_all_mr_jackknifed_INS_props(*jp_p.first,*jp_p.second,ijack);
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
  
  clusterize_all_mr_INS_props(jprop_0,clust_size);
  if(use_QED)
    for(auto &p : {&jprop_2,&jprop_P,&jprop_S})
      clusterize_all_mr_INS_props(*p,clust_size);
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

