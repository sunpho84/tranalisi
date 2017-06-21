#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <Zbil.hpp>

#include <prop.hpp>

void set_mr_Zbil_ind(size_t nm,size_t nr)
{
  mr_Zbil_ind.set_ranges({{"mr_bw",nmr},{"mr_fw",nmr},{"Zbil",nZbil}});
  nmr_Zbil=mr_Zbil_ind.max();
}

void set_mr_gbil_ind(size_t nm,size_t nr)
{
  mr_gbil_ind.set_ranges({{"mr_bw",nmr},{"mr_fw",nmr},{"iG",nGamma}});
  nmr_gbil=mr_gbil_ind.max();
}

void set_jbil_verts(bool set_QED)
{
  jverts_0.resize(nmr_gbil,vjbil_vert_t(imoms.size()));
  if(set_QED)
    for(auto &o : {&jverts_em,&jverts_P,&jverts_S})
      o->resize(nmr_gbil,vjbil_vert_t(imoms.size()));
}

void build_all_mr_gbil_jackknifed_INS_verts(vector<vjbil_vert_t> &jverts,const vector<vprop_t> &prop1,const vector<vprop_t> &prop2,size_t ijack)
{
  for(size_t imr_bw=0;imr_bw<nmr;imr_bw++)
    for(size_t imr_fw=0;imr_fw<nmr;imr_fw++)
      for(size_t imom=0;imom<imoms.size();imom++)
	for(size_t iG=0;iG<nGamma;iG++)
	  add_to_cluster(jverts[mr_gbil_ind({imr_bw,imr_fw,iG})][imom],prop1[imr_fw][imom]*Gamma[iG]*Gamma[5]*prop2[imr_bw][imom].adjoint()*Gamma[5],ijack);
}

void build_all_mr_gbil_jackknifed_verts(bool use_QED,size_t ijack)
{
  build_all_mr_gbil_jackknifed_INS_verts(jverts_0,conf_prop_0,conf_prop_0,ijack);
  
  if(use_QED)
  for(auto &o : vector<tuple<vector<vjbil_vert_t>*,vector<vprop_t>*,vector<vprop_t>*>>{
    {&jverts_em,&conf_prop_F,&conf_prop_F},
    {&jverts_em,&conf_prop_FF,&conf_prop_0},
    {&jverts_em,&conf_prop_0,&conf_prop_FF},
    {&jverts_em,&conf_prop_T,&conf_prop_0},
    {&jverts_em,&conf_prop_0,&conf_prop_T},
    {&jverts_P,&conf_prop_P,&conf_prop_0},
    {&jverts_P,&conf_prop_0,&conf_prop_P},
    {&jverts_S,&conf_prop_S,&conf_prop_0},
    {&jverts_S,&conf_prop_0,&conf_prop_S}})
    build_all_mr_gbil_jackknifed_INS_verts(*get<0>(o),*get<1>(o),*get<2>(o),ijack);
}

void clusterize_all_mr_gbil_INS_verts(vector<vjbil_vert_t> &jverts,size_t clust_size)
{
  for(size_t imr_gbil=0;imr_gbil<nmr_gbil;imr_gbil++)
#pragma omp parallel for
    for(size_t imom=0;imom<imoms.size();imom++)
      clusterize(jverts[imr_gbil][imom],clust_size);
}

void clusterize_all_mr_gbil_verts(bool use_QED,size_t clust_size)
{
  clusterize_all_mr_gbil_INS_verts(jverts_0,clust_size);
  if(use_QED)
    for(auto &v : {&jverts_em,&jverts_P,&jverts_S})
      clusterize_all_mr_gbil_INS_verts(*v,clust_size);
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

