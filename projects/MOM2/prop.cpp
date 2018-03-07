#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PROP
 #include <MOM2/prop.hpp>

#include <MOM2/pars.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

string get_qprop_tag(const size_t im,const size_t ir,const size_t ikind)
{
  return combine("S_M%zu_R%zu_%s",im,ir,m_r_mom_conf_qprops_t::tag[ikind]);
}

string get_lprop_tag(const size_t ikind)
{
  return combine("L_%s",mom_conf_lprops_t::tag[ikind]);
}

vector<raw_file_t> perens_t::setup_read_all_qprops_mom(const vector<size_t> &conf_list) const
{
  using namespace pars;
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
  vector<raw_file_t> files(im_r_iconf_ihit_iqkind_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<im_r_iconf_ihit_iqkind_ind.max();i++)
    {
      const vector<size_t> comps=im_r_iconf_ihit_iqkind_ind(i);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t iconf=comps[2];
      const size_t ihit=comps[3];
      const size_t ikind=comps[4];
      const string path_base=combine("%s/%04zu/fft_",prop_hadr_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      
      files[i].open(path_base+get_qprop_tag(im,r,ikind)+path_suff,"r");
    }
  
  return files;
}

vector<raw_file_t> perens_t::setup_read_all_lprops_mom(const vector<size_t> &conf_list) const
{
  vector<raw_file_t> files(iconf_ihit_ilkind_ind.max());
  
  string suff_hit="";
  if(nhits>1) suff_hit="_hit_%zu";
  
#pragma omp parallel for
  for(size_t i=0;i<iconf_ihit_ilkind_ind.max();i++)
    {
      const vector<size_t> comps=iconf_ihit_ilkind_ind(i);
      const size_t iconf=comps[0];
      const size_t ihit=comps[1];
      const size_t ikind=comps[2];
      const string path_base=combine("%s/%04zu/fft_",prop_lep_path.c_str(),conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      
      files[i].open(path_base+get_lprop_tag(ikind)+path_suff,"r");
    }
  
  return files;
}


void perens_t::mom_compute_qprop()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
    {
      const size_t mom=linmoms[ilinmom][0];
      vector<jm_r_mom_qprops_t> jprops(im_r_ind.max()); //!< jackknived props
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    const size_t mom=linmoms[ilinmom][0];
	    cout<<"Working on qprop, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum "<<ilinmom+1<<"/"<<linmoms.size()<<", "
	      "mom: "<<mom<<endl;
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props=read_all_qprops_mom(files,im_r_ijack_ind,i_in_clust_ihit,mom);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops,props,pars::use_QED,im_r_ind);
	    build_props_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops,pars::use_QED,clust_size,im_r_ind,deltam_cr);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv; //!< inverse propagator
      vector<jqprop_t> jprop_QED_inv; //!< inverse propagator with em insertion
      
      get_inverse_propagators(jprop_inv,jprop_QED_inv,jprops,im_r_ijackp1_ind);
      
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=im_r_ind({im,r});
	  const size_t im_r_ilinmom=im_r_ilinmom_ind({im,r,ilinmom});
	  
	  //compute Zq
	  Zq_time.start();
	  Zq[im_r_ilinmom][ijack]=compute_Zq(jprop_inv[im_r][ijack],mom);
	  Zq_sig1[im_r_ilinmom][ijack]=compute_Zq_sig1(jprop_inv[im_r][ijack],mom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(pars::use_QED)
	    {
	      Zq_time.start();
	      Zq_QED[im_r_ilinmom][ijack]=-compute_Zq(jprop_QED_inv[im_r][ijack],mom);
	      Zq_sig1_QED[im_r_ilinmom][ijack]=-compute_Zq_sig1(jprop_QED_inv[im_r][ijack],mom);
	      Zq_time.stop();
	    }
	}
    }
}

void perens_t::get_inverse_propagators(vector<jqprop_t> &jprop_inv,vector<jqprop_t> &jprop_QED_inv,
				       const vector<jm_r_mom_qprops_t> &jprops,
				       const index_t &im_r_ijackp1_ind) const
{
  jprop_inv.resize(im_r_ind.max());
  jprop_QED_inv.resize(im_r_ind.max());
  
#pragma omp parallel for reduction(+:invert_time)
  for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
    {
      //decript indices
      const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
      const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
      const size_t im_r=im_r_ind({im,r});
      
      //compute inverse
      invert_time.start();
      qprop_t prop_inv=jprop_inv[im_r][ijack]=jprops[im_r].LO[ijack].inverse();
      invert_time.stop();
      
      //do the same with QED
      if(pars::use_QED)
	{
	  invert_time.start(); //This misses a sign -1 coming from the original inverse
	  jprop_QED_inv[im_r][ijack]=prop_inv*jprops[im_r].QED[ijack]*prop_inv;
	  invert_time.stop();
	}
    }
}

