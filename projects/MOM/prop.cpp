#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

#include <prop.hpp>

vjprop_t get_prop_inv(const vjprop_t &jprop)
{
  vjprop_t jprop_inv(jprop.size());
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      put_into_jackknife(jprop_inv[imom],get_from_jackknife(jprop[imom],ijack).inverse(),ijack);
  
  return jprop_inv;
}

vprop_t read_prop(const string &template_path,size_t iconf)
{
  vprop_t prop(imoms.size());
  
  //! source file
  raw_file_t file(combine(template_path.c_str(),iconf),"r");
  
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t imom=0;imom<imoms.size();imom++)
	for(size_t is_si=0;is_si<NSPIN;is_si++)
	  for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	    file.bin_read(prop[imom](isc(is_si,ic_si),isc(is_so,ic_so)));
  
  return prop;
}

void build_jackknifed_prop(vjprop_t &jprop,const vprop_t &prop,size_t ijack)
{
  for(size_t imom=0;imom<imoms.size();imom++)
    add_to_cluster(jprop[imom],prop[imom],ijack);
}

void clusterize_prop(vjprop_t &jprop,size_t clust_size)
{
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    clusterize(jprop[imom],clust_size);
}
