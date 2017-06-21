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

vprop_t read_prop(const string &path)
{
  vprop_t prop(imoms.size());
  
  //! source file
  raw_file_t file(path,"r");
  
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t iall_mom=0,imom=0;iall_mom<filt_moms.size();iall_mom++)
	{
	  for(size_t is_si=0;is_si<NSPIN;is_si++)
	    for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	      {
		dcompl_t c;
		file.bin_read(c);
		if(filt_moms[iall_mom]) prop[imom](isc(is_si,ic_si),isc(is_so,ic_so))=c;
	    }
	  if(filt_moms[iall_mom]) imom++;
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
  conf_prop_0.resize(nmr);
  if(set_QED)
    {
      conf_prop_FF.resize(nmr);
      conf_prop_F.resize(nmr);
      conf_prop_T.resize(nmr);
      conf_prop_P.resize(nmr);
      conf_prop_S.resize(nmr);
    }
}

void read_all_mr_INS_props(vector<vprop_t> &conf_prop,const string &template_path,const string &ins)
{
  for(size_t im=0;im<nm;im++)
    for(size_t r=0;r<nr;r++)
      conf_prop[mr_ind({im,r})]=read_prop(combine(template_path.c_str(),get_prop_tag(im,r,ins).c_str()));
}

void read_all_mr_QED_props(const string &template_path)
{
  for(auto &o : vector<pair<vector<vprop_t>*,string>>({{&conf_prop_FF,"FF"},{&conf_prop_F,"F"},{&conf_prop_T,"T"},{&conf_prop_P,"P"},{&conf_prop_S,""}}))
    read_all_mr_INS_props(*o.first,template_path,o.second);
  
  dcompl_t fact_P(0.0,-1.0);
  dcompl_t fact_S(-1.0,0.0);
  
  //put factors
  for(auto &pmr_f : vector<pair<vector<vprop_t>*,dcompl_t>>({{&conf_prop_P,fact_P},{&conf_prop_S,fact_S}}))
    for(auto &p : *pmr_f.first) //all mr
      for(auto &p_mom : p) //all momentum
	p_mom*=pmr_f.second;
}

void read_all_mr_props(bool read_QED,const string &template_path)
{
  read_all_mr_LO_props(template_path);
  if(read_QED) read_all_mr_QED_props(template_path);
}

void set_jprops(bool set_QED)
{
  cout<<"Setting all "<<nmr<<" jprops"<<endl;
  
  jprop_0.resize(nmr,vjprop_t(imoms.size()));
  if(set_QED)
    for(auto &o : {&jprop_2,&jprop_P,&jprop_S})
      o->resize(nmr,vjprop_t(imoms.size()));
}

void build_all_mr_jackknifed_INS_props(vector<vjprop_t> &out,const vector<vprop_t> &in,size_t ijack)
{
  for(size_t imr=0;imr<nmr;imr++)
    for(size_t imom=0;imom<imoms.size();imom++)
      add_to_cluster(out[imr][imom],in[imr][imom],ijack);
}

void build_all_mr_jackknifed_props(bool set_QED,size_t ijack)
{
  build_all_mr_jackknifed_INS_props(jprop_0,conf_prop_0,ijack);
  if(set_QED)
    for(auto &jp_p : vector<pair<vector<vjprop_t>*,vector<vprop_t>*>>({{&jprop_2,&conf_prop_FF},{&jprop_2,&conf_prop_T},{&jprop_P,&conf_prop_P},{&jprop_S,&conf_prop_S}}))
  build_all_mr_jackknifed_INS_props(*jp_p.first,*jp_p.second,ijack);
}

void clusterize_all_mr_INS_props(vector<vjprop_t> &jprop,size_t clust_size)
{
  for(size_t imr=0;imr<nmr;imr++)
#pragma omp parallel for
    for(size_t imom=0;imom<imoms.size();imom++)
      clusterize(jprop[imr][imom],clust_size);
}

void clusterize_all_mr_props(bool use_QED,size_t clust_size)
{
  clusterize_all_mr_INS_props(jprop_0,clust_size);
  if(use_QED)
    for(auto &p : {&jprop_2,&jprop_P,&jprop_S})
      clusterize_all_mr_INS_props(*p,clust_size);
}

vector<vjprop_t> get_all_mr_props_inv(const vector<vjprop_t> &jprop)
{
  vector<vjprop_t> jprop_inv(jprop.size());
  
  for(size_t imr=0;imr<nmr;imr++)
#pragma omp parallel for
    for(size_t imom=0;imom<imoms.size();imom++)
      for(size_t ijack=0;ijack<=njacks;ijack++)
	put_into_jackknife(jprop_inv[imr][imom],get_from_jackknife(jprop[imr][imom],ijack).inverse(),ijack);
  
  return jprop_inv;
}

