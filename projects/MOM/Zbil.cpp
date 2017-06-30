#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <Zbil.hpp>

#include <prop.hpp>

void build_jackknifed_vert_Gamma(jprop_t &jvert,const prop_t &prop1,size_t iG,const prop_t &prop2,size_t ijack)
{add_to_cluster(jvert,prop1*Gamma[iG]*Gamma[5]*prop2.adjoint()*Gamma[5],ijack);}

void build_all_mr_gbil_jackknifed_verts(vector<jbil_vert_t> &jbil,const vector<m_r_mom_conf_props_t> &props,
					const index_t &im_r_im_r_igam_ind,const index_t &im_r_ijack_ind,bool use_QED)
{
  //! help finding the bilinear/jack combo
  index_t ind({{"i",im_r_im_r_igam_ind.max()},{"ijack",njacks}});
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //decript bilinar/jack
      vector<size_t> ind_comp=ind(i);
      const size_t im_r_im_r_igam=ind_comp[0],ijack=ind_comp[1];
      
      //decript props to combine and which Gamma to create
      const vector<size_t> im_r_im_r_igam_comp=im_r_im_r_igam_ind(im_r_im_r_igam);
      const size_t im_fw=im_r_im_r_igam_comp[0],r_fw=im_r_im_r_igam_comp[1];
      const size_t im_bw=im_r_im_r_igam_comp[2],r_bw=im_r_im_r_igam_comp[3];
      const size_t iG=im_r_im_r_igam_comp[4];
      
      //proxy for vector and props
      jbil_vert_t &jv=jbil[im_r_im_r_igam];
      const m_r_mom_conf_props_t &p1=props[im_r_ijack_ind({im_fw,r_fw,ijack})];
      const m_r_mom_conf_props_t &p2=props[im_r_ijack_ind({im_bw,r_bw,ijack})];
      
      //create list of operations
      vector<tuple<jprop_t*,const prop_t*,const prop_t*>> list={{&jv.jvert_0,&p1.prop_0,&p2.prop_0}};
      if(use_QED)
	for(auto &o : vector<tuple<jprop_t*,const prop_t*,const prop_t*>>({
	    {&jv.jvert_em,&p1.prop_F,&p2.prop_F},
	    {&jv.jvert_em,&p1.prop_FF,&p2.prop_0},
            {&jv.jvert_em,&p1.prop_0,&p2.prop_FF},
            {&jv.jvert_em,&p1.prop_T,&p2.prop_0},
            {&jv.jvert_em,&p1.prop_0,&p2.prop_T},
            {&jv.jvert_P,&p1.prop_P,&p2.prop_0},
            {&jv.jvert_P,&p1.prop_0,&p2.prop_P},
            {&jv.jvert_S,&p1.prop_S,&p2.prop_0},
            {&jv.jvert_S,&p1.prop_0,&p2.prop_S}}))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackknifed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),ijack);
    }
}

void clusterize_all_mr_gbil_verts(vector<jbil_vert_t> &jbils,bool use_QED,size_t clust_size)
{
  cout<<"Clusterizing all gbil verts, clust_size="<<clust_size<<endl;
  
#pragma omp parallel for
  for(size_t ibil=0;ibil<jbils.size();ibil++)
    jbils[ibil].clusterize_all(use_QED,clust_size);
}

// vector<djvec_t> compute_proj_bil(const vjprop_t &jprop_inv1,vjbil_vert_t &jverts,const vjprop_t &jprop_inv2)
// {
//   const size_t iZbil_of_iG[nGamma]={iZS,iZA,iZA,iZA,iZA,iZP,iZV,iZV,iZV,iZV,iZT,iZT,iZT,iZT,iZT,iZT};
//   const double Zdeg[nZbil]={1,4,1,4,6};
//   vector<djvec_t> pr(nZbil,djvec_t(imoms.size(),0.0)); //!< Five summed projectors
  
// #pragma omp parallel for
//   for(size_t imom=0;imom<imoms.size();imom++)
//     for(size_t ijack=0;ijack<=njacks;ijack++)
//       {
// 	prop_t prop_inv1=get_from_jackknife(jprop_inv1[imom],ijack);
// 	prop_t prop_inv2=get_from_jackknife(jprop_inv2[imom],ijack);
	
// 	for(size_t iG=0;iG<nGamma;iG++)
// 	  {
// 	    size_t iZbil=iZbil_of_iG[iG];
// 	    prop_t vert=get_from_jackknife(jverts[imom][iG],ijack) ;
// 	    prop_t amp_vert=prop_inv1*vert*Gamma[5]*prop_inv2.adjoint()*Gamma[5];
// 	    pr[iZbil][imom][ijack]+=(amp_vert*Gamma[iG].adjoint()).trace().real()/(12.0*Zdeg[iZbil]);
// 	  }
//       }
  
//   return pr;
// }

