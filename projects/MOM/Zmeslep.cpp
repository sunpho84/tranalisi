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

void build_all_mr_gmeslep_jackknifed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
					   const vector<mom_conf_lprops_t> &props_lep,
					   const index_t &im_r_im_r_igam_iprojlep_ind,const index_t &im_r_iclust_ind,const double q1,const double q2)
{
  //! help finding the meslep/clust combo
  index_t ind({{"i",im_r_im_r_igam_iprojlep_ind.max()},{"iclust",njacks}});
  
  //! takes the correct element in the lepton insertion-projection
  index_t iG_iprojlep_iclust_ind({{"iG",nGamma},{"iprojlep",nGamma},{"iclust",njacks}});
  
  //! projected lepton propagator with insertion
  vector<dcompl_t> mesloop(iG_iprojlep_iclust_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<iG_iprojlep_iclust_ind.max();i++)
    {
      vector<size_t> ind_comp=iG_iprojlep_iclust_ind(i);
      const size_t iG=ind_comp[0],iprojlep=ind_comp[1],iclust=ind_comp[2];
      const mom_conf_lprops_t &pl=props_lep[iclust];
      mesloop[i]=(pl.F*lepGamma[iG]*(lepGamma[0]-lepGamma[5])*lepGamma[iprojlep]).trace()/4.0;
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
	 {{&j.ML1[im_r_im_r_igam_iprojlep],&p1.F,&p2.LO,q1*mesloop[imesloop]},
	  {&j.ML2[im_r_im_r_igam_iprojlep],&p1.LO,&p2.F,q2*mesloop[imesloop]}};
      
      //create the vertex
      for(auto &o : list)
	build_jackknifed_meslep_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),get<3>(o),iclust);
    }
}

// djvec_t compute_proj_bil(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind)
// {
//   const size_t nm=im_r_ind.max(0),nr=im_r_ind.max(1);
//   const index_t im_r_im_r_iG_ind({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma",nGamma}});
//   const index_t im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"Zbil",nZbil}});
//   const index_t ind({{"rest",im_r_im_r_iZbil_ind.max()},{"ijack",njacks+1}});
  
//   //combine the 16 bilinears to form the 5 vertex
//   const vector<vector<size_t>> iG_of_Zbil={{0},{1,2,3,4},{5},{6,7,8,9},{10,11,12,13,14,15}};
//   djvec_t pr(im_r_im_r_iZbil_ind.max());
  
// #pragma omp parallel for
//   for(size_t i=0;i<ind.max();i++)
//     {
//       //split im_r_im_r_iZbil and ijack
//       vector<size_t> i_comp=ind(i);
//       const size_t im_r_im_r_iZbil=i_comp[0],ijack=i_comp[1];
//       const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
      
//       //get im and r for fw and back, and iZbil
//       const size_t im_fw=im_r_im_r_iZbil_comp[0],r_fw=im_r_im_r_iZbil_comp[1];
//       const size_t im_bw=im_r_im_r_iZbil_comp[2],r_bw=im_r_im_r_iZbil_comp[3];
//       const size_t iZbil=im_r_im_r_iZbil_comp[4];
      
//       //loop on all gammas
//       djack_t &out=pr[im_r_im_r_iZbil];
//       for(auto & iG : iG_of_Zbil[iZbil])
// 	{
// 	  vector<size_t> im_r_im_r_iG_comp=im_r_im_r_iZbil_comp;
// 	  im_r_im_r_iG_comp[4]=iG;
// 	  const size_t im_r_im_r_iG=im_r_im_r_iG_ind(im_r_im_r_iG_comp);
	  
// 	  const size_t ip1=im_r_ind({im_fw,r_fw});
// 	  const size_t ip2=im_r_ind({im_bw,r_bw});
// 	  const qprop_t &prop_inv1=jprop_inv1[ip1][ijack];
// 	  const qprop_t &prop_inv2=jprop_inv2[ip2][ijack];
	  
// 	  const qprop_t &vert=jverts[im_r_im_r_iG][ijack];
	  
// 	  qprop_t amp_vert=prop_inv1*vert*Gamma[5]*prop_inv2.adjoint()*Gamma[5];
// 	  out[ijack]+=(amp_vert*Gamma[iG].adjoint()).trace().real()/(12.0*iG_of_Zbil[iZbil].size());
// 	}
//     }
  
//   return pr;
//
// }
