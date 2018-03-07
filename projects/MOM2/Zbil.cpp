#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <MOM2/Zbil.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

void perens_t::mom_compute_bil()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    {
      const size_t imom1=bilmoms[ibilmom][1];
      const size_t imom2=bilmoms[ibilmom][2];
      const bool read2=(imom1!=imom2);
      
      vector<jm_r_mom_qprops_t> jprops1(im_r_ind.max()); //!< jackknived props
      vector<jm_r_mom_qprops_t> jprops2(im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),pars::use_QED);   //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t mom1=linmoms[imom1][0];
	    const size_t mom2=linmoms[imom2][0];
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on qbil, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum combo "<<ibilmom+1<<"/"<<bilmoms.size()<<", "
	      "moms: "<<mom1<<" "<<mom2<<endl;
	    
	    //read
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props1=read_all_qprops_mom(files,im_r_ijack_ind,i_in_clust_ihit,mom1);
	    const vector<m_r_mom_conf_qprops_t> props2=(read2?read_all_qprops_mom(files,im_r_ijack_ind,i_in_clust_ihit,mom2):props1);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops1,props1,pars::use_QED,im_r_ind);
	    build_all_mr_jackkniffed_qprops(jprops2,props2,pars::use_QED,im_r_ind);
	    build_props_time.stop();
	    
	    //build all bilinear verts
	    build_verts_time.start();
	    build_all_mr_gbil_jackkniffed_verts(jverts,props1,props2,im_r_im_r_igam_ind,im_r_ijack_ind,pars::use_QED);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops1,pars::use_QED,clust_size,im_r_ind,deltam_cr);
      clusterize_all_mr_jackkniffed_qprops(jprops2,pars::use_QED,clust_size,im_r_ind,deltam_cr);
      jverts.clusterize_all(pars::use_QED,clust_size,im_r_im_r_igam_ind,deltam_cr);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv1; //!< inverse propagator1
      vector<jqprop_t> jprop_inv2; //!< inverse propagator2
      vector<jqprop_t> jprop_QED_inv1; //!< inverse propagator1 with em insertion
      vector<jqprop_t> jprop_QED_inv2; //!< inverse propagator2 with em insertion
      
      get_inverse_propagators(jprop_inv1,jprop_QED_inv1,jprops1,im_r_ijackp1_ind);
      get_inverse_propagators(jprop_inv2,jprop_QED_inv2,jprops2,im_r_ijackp1_ind);
      
      proj_time.start();
      djvec_t pr_bil_temp=compute_proj_bil(jprop_inv1,jverts.LO,jprop_inv2,im_r_ind);
      
      //QED
      djvec_t pr_bil_QED_temp;
      if(pars::use_QED)
	{
	  const djvec_t pr_bil_QED=compute_proj_bil(jprop_inv1,jverts.QED,jprop_inv2,im_r_ind);
	  const djvec_t pr_bil_a=compute_proj_bil(jprop_QED_inv1,jverts.LO,jprop_inv2,im_r_ind);
	  const djvec_t pr_bil_b=compute_proj_bil(jprop_inv1,jverts.LO,jprop_QED_inv2,im_r_ind);
	  pr_bil_QED_temp=
	    -pr_bil_a
	    -pr_bil_b
	    +pr_bil_QED;
	}
      
      //! an index running on all packed combo, and momenta
      const index_t all_ibilmom_ind({{"All",im_r_im_r_iZbil_ind.max()},{"bilmom",bilmoms.size()}});
      
      //store
      vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil,&pr_bil_temp}};
      if(pars::use_QED) store_task.push_back({&pr_bil_QED,&pr_bil_QED_temp});
      for(auto &t : store_task)
	for(size_t iall=0;iall<im_r_im_r_iZbil_ind.max();iall++)
	  (*t.first)[all_ibilmom_ind({iall,ibilmom})]=(*t.second)[iall];
      
      proj_time.stop();
    }
}
