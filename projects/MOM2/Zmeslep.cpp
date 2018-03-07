#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_MESLEP
 #include <MOM2/Zmeslep.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

void perens_t::mom_compute_meslep()
{
  vector<raw_file_t> qfiles=setup_read_all_qprops_mom(conf_list);
  vector<raw_file_t> lfiles=setup_read_all_lprops_mom(conf_list);
  
  for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
    {
      const size_t imom_in=meslepmoms()[imeslepmom][1];
      const size_t imom_ou=meslepmoms()[imeslepmom][2];
      const bool read_ou=(imom_in!=imom_ou);
      
      vector<jm_r_mom_qprops_t> jprops_in(im_r_ind.max());                //!< jackknived props
      vector<jm_r_mom_qprops_t> jprops_ou(im_r_ind.max());                //!< jackknived props
      jmeslep_vert_t jmeslep_verts(im_r_im_r_ilistGl_ipGl_ind.max());     //!< jackknived vertex
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t mom_in=linmoms[imom_in][0];
	    const size_t mom_ou=linmoms[imom_ou][0];
	    const size_t momlep=mom_in;
	    if(pars::scheme==reno_scheme::SMOM) CRASH("this will not work for SMOM");
	    
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on meslep, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum combo "<<imeslepmom+1<<"/"<<meslepmoms().size()<<", "
	      "moms: "<<mom_in<<" "<<mom_ou<<endl;
	    
	    //read
	    read_time.start();
	    vector<m_r_mom_conf_qprops_t> props_in=read_all_qprops_mom(qfiles,im_r_ijack_ind,i_in_clust_ihit,mom_in);
	    vector<m_r_mom_conf_qprops_t> props_ou=(read_ou?read_all_qprops_mom(qfiles,im_r_ijack_ind,i_in_clust_ihit,mom_ou):props_in);
	    vector<mom_conf_lprops_t> props_lep=read_all_lprops_mom(lfiles,i_in_clust_ihit,momlep);
	    read_time.stop();
	    
	    //put the charges in place
	    using namespace meslep;
	    incorporate_charge(props_in,q_in);
	    incorporate_charge(props_ou,q_ou);
	    incorporate_charge(props_lep,ql);
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops_in,props_in,pars::use_QED,im_r_ind);
	    build_all_mr_jackkniffed_qprops(jprops_ou,props_ou,pars::use_QED,im_r_ind);
	    build_props_time.stop();
	    
	    //build all meslep_verts
	    build_meslep_verts_time.start();
	    build_all_mr_gmeslep_jackkniffed_verts(jmeslep_verts,props_in,props_ou,props_lep,im_r_ind);
	    build_meslep_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops_in,pars::use_QED,clust_size,im_r_ind,deltam_cr);
      clusterize_all_mr_jackkniffed_qprops(jprops_ou,pars::use_QED,clust_size,im_r_ind,deltam_cr);
      jmeslep_verts.clusterize_all(clust_size,im_r_im_r_iop_ilistpGl_ind,deltam_cr);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv_in;      //!< inverse propagator_in
      vector<jqprop_t> jprop_inv_ou;      //!< inverse propagator_ou
      vector<jqprop_t> jprop_QED_inv_in;  //!< inverse propagator_in with em insertion
      vector<jqprop_t> jprop_QED_inv_ou;  //!< inverse propagator_ou with em insertion
      
      get_inverse_propagators(jprop_inv_in,jprop_QED_inv_in,jprops_in,im_r_ijackp1_ind);
      get_inverse_propagators(jprop_inv_ou,jprop_QED_inv_ou,jprops_ou,im_r_ijackp1_ind);
      
      proj_time.start();
      
      const auto &j=jmeslep_verts;
      
      djvec_t pr_LO,pr_QED_amp_QCD,pr_QCD_amp_QED_in,pr_QCD_amp_QED_ou;
      for(auto &p : vector<tuple<djvec_t*,const vector<jqprop_t>*,const vector<jqprop_t>*,const vector<jqprop_t>*>>{
	  {&pr_LO,              &jprop_inv_in,     &j.LO,  &jprop_inv_ou},
	  {&pr_QED_amp_QCD,     &jprop_inv_in,     &j.QED, &jprop_inv_ou},
	  {&pr_QCD_amp_QED_in,  &jprop_QED_inv_in, &j.LO,  &jprop_inv_ou},
	  {&pr_QCD_amp_QED_ou,  &jprop_inv_in,     &j.LO,  &jprop_QED_inv_ou}
	    })
	{
	  auto &out=*get<0>(p);
	  auto &p_in=*get<1>(p);
	  auto &v=*get<2>(p);
	  auto &p_ou=*get<3>(p);
	  out=compute_proj_meslep(p_in,v,p_ou,im_r_ind);
	}
      proj_time.stop();
      
      //! an index running on all packed combo, and momenta
      const index_t all_imeslepmom_ind({{"All",im_r_im_r_iop_iproj_ind.max()},{"meslepmom",meslepmoms().size()}});
      
      //store
      for(size_t iall=0;iall<im_r_im_r_iop_iproj_ind.max();iall++)
	{
	  pr_meslep[all_imeslepmom_ind({iall,imeslepmom})]=
	    +pr_LO[iall];
	  pr_meslep_QED[all_imeslepmom_ind({iall,imeslepmom})]=
	    +pr_QED_amp_QCD[iall]
	    -pr_QCD_amp_QED_in[iall]
	    -pr_QCD_amp_QED_ou[iall];
	}
    }
}
