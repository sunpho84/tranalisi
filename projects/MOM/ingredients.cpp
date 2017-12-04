#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_INGREDIENTS
 #include <ingredients.hpp>

#include <deltam_cr.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

//! Constructor
ingredients_t::ingredients_t()
{
  Zq_allmoms.resize(im_r_imom_ind.max());
  Zq_sig1_allmoms.resize(im_r_imom_ind.max());
  if(use_QED) Zq_sig1_EM_allmoms.resize(im_r_imom_ind.max());
  
  pr_bil_mom.resize(im_r_im_r_iZbil_imom_ind.max());
  if(use_QED) pr_bil_QED_mom.resize(im_r_im_r_iZbil_imom_ind.max());
}

void ingredients_t::ri_mom()
{
  vector<m_r_mom_conf_props_t> props; //!< store props for individual conf
  
  vector<task_list_t> read_tasks=prepare_read_prop_taks(props,conf_list);
  
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      vector<jm_r_mom_props_t> jprops(im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED); //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_hit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", hit "<<ihit+1<<"/"<<nhits<<", momentum "<<imom+1<<"/"<<imoms.size()<<endl;
	    read_time.start();
	    read_tasks[i_in_clust_hit].assolve_all(RECYCLE);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackknifed_props(jprops,props,use_QED,im_r_ind,deltam_cr);
	    build_props_time.stop();
	    
	    build_verts_time.start();
	    build_all_mr_gbil_jackknifed_verts(jverts,props,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED,deltam_cr);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
      clust_time.stop();
      
      vector<jprop_t> jprop_inv(im_r_ind.max()); //!< inverse propagator
      vector<jprop_t> jprop_EM_inv(im_r_ind.max()); //!< inverse propagator with em insertion
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=im_r_ind({im,r});
	  const size_t im_r_imom=im_r_imom_ind({im,r,imom});
	  
	  //compute inverse
	  invert_time.start();
	  prop_t prop_inv=jprops[im_r].LO[ijack].inverse();
	  //if(im_r_ijack==0) cout<<jprops[im_r].LO[ijack](0,0)<<endl;
	  jprop_inv[im_r][ijack]=prop_inv;
	  invert_time.stop();
	  
	  //compute Zq
	  Zq_time.start();
	  Zq_allmoms[im_r_imom][ijack]=compute_Zq(prop_inv,imom);
	  Zq_sig1_allmoms[im_r_imom][ijack]=compute_Zq_sig1(prop_inv,imom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(use_QED)
	    {
	      invert_time.start();
	      prop_t prop_EM_inv=prop_inv*jprops[im_r].EM[ijack]*prop_inv;
	      jprop_EM_inv[im_r][ijack]=prop_EM_inv;
	      invert_time.stop();
	      
	      Zq_time.start();
	      Zq_sig1_EM_allmoms[im_r_imom][ijack]=-compute_Zq_sig1(prop_EM_inv,imom);
	      Zq_time.stop();
	    }
	}
      
      //! an index running on all packed combo, and momenta
      index_t all_imom_ind({{"All",im_r_im_r_iZbil_ind.max()},{"mom",imoms.size()}});
      
      proj_time.start();
      djvec_t pr_bil_mom_temp=compute_proj_bil(jprop_inv,jverts.LO,jprop_inv,im_r_ind);
      
      //QED
      djvec_t pr_bil_QED_mom_temp;
      if(use_QED)
	{
	  djvec_t pr_bil_EM_mom=compute_proj_bil(jprop_inv,jverts.EM,jprop_inv,im_r_ind);
	  djvec_t pr_bil_a_mom=compute_proj_bil(jprop_EM_inv,jverts.LO,jprop_inv,im_r_ind);
	  djvec_t pr_bil_b_mom=compute_proj_bil(jprop_inv,jverts.LO,jprop_EM_inv,im_r_ind);
	  pr_bil_QED_mom_temp=pr_bil_a_mom+pr_bil_b_mom-pr_bil_EM_mom;
	}
      
      //store
      vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil_mom,&pr_bil_mom_temp}};
      if(use_QED) store_task.push_back({&pr_bil_QED_mom,&pr_bil_QED_mom_temp});
      for(auto &t : store_task)
	for(size_t iall=0;iall<im_r_im_r_iZbil_ind.max();iall++)
	  (*t.first)[all_imom_ind({iall,imom})]=(*t.second)[iall];
      
      proj_time.stop();
    }
}

void ingredients_t::bin_read(raw_file_t &file)
{
  for(djvec_t *o: {&Zq_allmoms,&Zq_sig1_allmoms,&Zq_sig1_EM_allmoms,&pr_bil_mom,&pr_bil_QED_mom}) o->bin_read(file);
}

void ingredients_t::bin_write(raw_file_t &file) const
{
  for(const djvec_t *o: {&Zq_allmoms,&Zq_sig1_allmoms,&Zq_sig1_EM_allmoms,&pr_bil_mom,&pr_bil_QED_mom}) o->bin_write(file);
}


