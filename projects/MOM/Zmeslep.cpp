#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_MESLEP
 #include <Zmeslep.hpp>

#include <ingredients.hpp>
#include <oper.hpp>

#include <prop.hpp>

void build_jackknifed_meslep_vert_Gamma(jqprop_t &jvert,const qprop_t &prop1,size_t iG,const qprop_t &prop2,const dcompl_t &lloop,size_t ijack)
{
  jvert[ijack]+=prop1*quaGamma[iG]*quaGamma[5]*prop2.adjoint()*quaGamma[5]*lloop;
}

void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					   const vector<mom_conf_lprops_t> &props_lep,
					   const index_t &im_r_im_r_igam_iprojlep_ind,const index_t &im_r_iclust_ind)
{
  //! help finding the meslep/clust combo
  index_t ind({{"i",im_r_im_r_igam_iprojlep_ind.max()},{"iclust",njacks}});
  
  //! takes the correct element in the lepton insertion-projection
  index_t iG_iprojlep_iclust_ind({{"iG",nGamma},{"iprojlep",nGamma},{"iclust",njacks}});
  
  //! projected lepton propagator with insertion
  vector<dcompl_t> mesloop(iG_iprojlep_iclust_ind.max());
  
  // NB: the lepton loop is fully amputated
#pragma omp parallel for
  for(size_t i=0;i<iG_iprojlep_iclust_ind.max();i++)
    {
      vector<size_t> ind_comp=iG_iprojlep_iclust_ind(i);
      const size_t iG=ind_comp[0],iprojlep=ind_comp[1],iclust=ind_comp[2];
      const mom_conf_lprops_t &pl=props_lep[iclust];
      mesloop[i]=(pl.F*lepGamma[iG]*(lepGamma[0]-lepGamma[5])*lepGamma[iprojlep]).trace()/4.0; //normalization for the single gamma
    }
  
 #pragma omp parallel for
   for(size_t i=0;i<ind.max();i++)
     {
       //decript bilinar/jack
       vector<size_t> ind_comp=ind(i);
       const size_t im_r_im_r_igam_iprojlep=ind_comp[0],iclust=ind_comp[1];
       
       //decript props to combine and which Gamma to create
       const vector<size_t> im_r_im_r_igam_iprojlep_comp=im_r_im_r_igam_iprojlep_ind(im_r_im_r_igam_iprojlep);
       const size_t im_fw=im_r_im_r_igam_iprojlep_comp[0],r_fw=im_r_im_r_igam_iprojlep_comp[1];
       const size_t im_bw=im_r_im_r_igam_iprojlep_comp[2],r_bw=im_r_im_r_igam_iprojlep_comp[3];
       const size_t iG=im_r_im_r_igam_iprojlep_comp[4];
       const size_t iprojlep=im_r_im_r_igam_iprojlep_comp[5];
       
       //proxy for vector and props
       const m_r_mom_conf_qprops_t &p1=props1[im_r_iclust_ind({im_fw,r_fw,iclust})];
       const m_r_mom_conf_qprops_t &p2=props2[im_r_iclust_ind({im_bw,r_bw,iclust})];
       
       //create list of operations
       const size_t imesloop=iG_iprojlep_iclust_ind({iG,iprojlep,iclust});
       vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*,const dcompl_t>> list=
	 {{&j.ML1[im_r_im_r_igam_iprojlep],&p1.F,&p2.LO,mesloop[imesloop]},
	  {&j.ML2[im_r_im_r_igam_iprojlep],&p1.LO,&p2.F,mesloop[imesloop]}};
      
      //create the vertex
      for(auto &o : list)
	build_jackknifed_meslep_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),get<3>(o),iclust);
    }
}

djvec_t compute_proj_measlep(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind)
{
  const size_t nm=im_r_ind.max(0),nr=im_r_ind.max(1);
  const index_t im_r_im_r_iGq_iPl_ind({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma_quark",nGamma},{"iproj_lep",nGamma}});
  const index_t im_r_im_r_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"iop",nZbil},{"iproj",nZbil}});
  const index_t ind({{"rest",im_r_im_r_iop_iproj_ind.max()},{"ijack",njacks+1}});
  
  //Each operator on the quark side has the structure G(1+-g5)
  //We decompose it in the 1 and g5 parts, putting explicitly the sign
  const vector<vector<size_t>> iGq_of_iop ={{ 1, 2, 3, 4}, { 1, 2, 3, 4}, { 0}, { 0}, {10,11,12,13,14,15}};
  const vector<int>            sign_of_iop={-1,             +1,            -1,   +1 ,  +1};
  const vector<vector<size_t>> &iPl_of_iproj=iGq_of_iop; // here instead 1-g5 was already included at lepton projection
  djvec_t pr(im_r_im_r_iop_iproj_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_iop_iproj and ijack
      vector<size_t> i_comp=ind(i);
      const size_t im_r_im_r_iop_iproj=i_comp[0],ijack=i_comp[1];
      const vector<size_t> im_r_im_r_iop_iproj_comp=im_r_im_r_iop_iproj_ind(im_r_im_r_iop_iproj);
      
      //get im and r for fw and back, and iop and iproj
      const size_t im_fw=im_r_im_r_iop_iproj_comp[0],r_fw=im_r_im_r_iop_iproj_comp[1];
      const size_t im_bw=im_r_im_r_iop_iproj_comp[2],r_bw=im_r_im_r_iop_iproj_comp[3];
      const size_t iop=im_r_im_r_iop_iproj_comp[4],iproj=im_r_im_r_iop_iproj_comp[5];
      
      //loop on all gammas
      djack_t &out=pr[im_r_im_r_iop_iproj];
      for(auto &iGq : iGq_of_iop[iop])
	{
 	  vector<size_t> im_r_im_r_iGq_iPl_comp=im_r_im_r_iop_iproj_comp;
	  im_r_im_r_iGq_iPl_comp[4]=iGq;
	  for(auto &iPl : iPl_of_iproj[iproj])
	    {
	      im_r_im_r_iGq_iPl_comp[5]=iPl;
	      const size_t im_r_im_r_iGq_iPl=im_r_im_r_iGq_iPl_ind(im_r_im_r_iGq_iPl_comp);
	      
	      const size_t ip1=im_r_ind({im_fw,r_fw});
	      const size_t ip2=im_r_ind({im_bw,r_bw});
	      const qprop_t &prop_inv1=jprop_inv1[ip1][ijack];
	      const qprop_t &prop_inv2=jprop_inv2[ip2][ijack];
	      
	      const qprop_t &vert=jverts[im_r_im_r_iGq_iPl][ijack];
	      
	      qprop_t amp_vert=prop_inv1*vert*quaGamma[5]*prop_inv2.adjoint()*quaGamma[5];
	      out[ijack]+=(amp_vert*quaGamma[iGq]*(quaGamma[0]+sign_of_iop[iop]*quaGamma[5]).adjoint()).trace().real()/(12.0*sqr(iGq_of_iop[iop].size()));
	    }
	}
    }
  
  return pr;
}
