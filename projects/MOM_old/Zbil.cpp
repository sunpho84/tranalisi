#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <Zbil.hpp>

#include <ingredients.hpp>
#include <oper.hpp>

#include <prop.hpp>

void build_jackkniffed_vert_Gamma(jqprop_t &jvert,const qprop_t &prop_in,size_t iG,const qprop_t &prop_ou,size_t iclust)
{
  auto c=prop_ou*quaGamma[iG]*quaGamma[5]*prop_in.adjoint()*quaGamma[5];
  jvert[iclust]+=c;
  
  if(0)
  cout
    <<"iG: "<<iG<<
    ", iclust: "<<iclust<<
    ", prop_in: "<<prop_in(0,0)<<
    ", prop_ou: "<<prop_ou(0,0)<<
    ", res: "<<c(0,0)<<
    endl;
}

void build_all_mr_gbil_jackkniffed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou,
					 const index_t &im_r_im_r_igam_ind,const index_t &im_r_iclust_ind,bool use_QED)
{
  //! help finding the bilinear/clust combo
  index_t ind({{"i",im_r_im_r_igam_ind.max()},{"iclust",njacks}});
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //decript bilinar/jack
      vector<size_t> ind_comp=ind(i);
      const size_t im_r_im_r_igam=ind_comp[0],iclust=ind_comp[1];
      
      //decript props to combine and which Gamma to create
      const vector<size_t> im_r_im_r_igam_comp=im_r_im_r_igam_ind(im_r_im_r_igam);
      const size_t im_fw=im_r_im_r_igam_comp[0],r_fw=im_r_im_r_igam_comp[1];
      const size_t im_bw=im_r_im_r_igam_comp[2],r_bw=im_r_im_r_igam_comp[3];
      const size_t iG=im_r_im_r_igam_comp[4];
      
      //proxy for vector and props
      const m_r_mom_conf_qprops_t &p_in=props_in[im_r_iclust_ind({im_fw,r_fw,iclust})];
      const m_r_mom_conf_qprops_t &p_ou=props_ou[im_r_iclust_ind({im_bw,r_bw,iclust})];
      
      //create list of operations
      vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>> list={{&jbil.LO[im_r_im_r_igam],&p_in.LO,&p_ou.LO}};
      if(use_QED)
	for(auto &o : vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>>({
	      {&jbil.PH[im_r_im_r_igam],&p_in.F,&p_ou.F},
// #warning including only exchange
	      {&jbil.PH[im_r_im_r_igam],&p_in.FF,&p_ou.LO},
	      {&jbil.PH[im_r_im_r_igam],&p_in.LO,&p_ou.FF},
	      {&jbil.PH[im_r_im_r_igam],&p_in.T,&p_ou.LO},
	      {&jbil.PH[im_r_im_r_igam],&p_in.LO,&p_ou.T},
	      {&jbil.CT_in[im_r_im_r_igam],&p_in.P,&p_ou.LO},
	      {&jbil.CT_ou[im_r_im_r_igam],&p_in.LO,&p_ou.P},
	      {&jbil.S[im_r_im_r_igam],&p_in.S,&p_ou.LO},
	      {&jbil.S[im_r_im_r_igam],&p_in.LO,&p_ou.S}
	    }))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackkniffed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),iclust);
    }
}

djvec_t compute_proj_bil(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou,const index_t &im_r_ind)
{
  const size_t nm=im_r_ind.max(0),nr=im_r_ind.max(1);
  const index_t im_r_im_r_iG_ind({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma",nGamma}});
  const index_t im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"Zbil",nZbil}});
  const index_t ind({{"rest",im_r_im_r_iZbil_ind.max()},{"ijack",njacks+1}});
  
  //combine the 16 bilinears to form the 5 vertex
  djvec_t pr(im_r_im_r_iZbil_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_iZbil and ijack
      vector<size_t> i_comp=ind(i);
      const size_t im_r_im_r_iZbil=i_comp[0],ijack=i_comp[1];
      const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
      
      //get im and r for fw and back, and iZbil
      const size_t im_fw=im_r_im_r_iZbil_comp[0],r_fw=im_r_im_r_iZbil_comp[1];
      const size_t im_bw=im_r_im_r_iZbil_comp[2],r_bw=im_r_im_r_iZbil_comp[3];
      const size_t iZbil=im_r_im_r_iZbil_comp[4];
      
      //loop on all gammas
      djack_t &out=pr[im_r_im_r_iZbil];
      for(auto & iG : iG_of_Zbil[iZbil])
	{
	  vector<size_t> im_r_im_r_iG_comp=im_r_im_r_iZbil_comp;
	  im_r_im_r_iG_comp[4]=iG;
	  const size_t im_r_im_r_iG=im_r_im_r_iG_ind(im_r_im_r_iG_comp);
	  
	  const size_t ip_in=im_r_ind({im_fw,r_fw});
	  const size_t ip_ou=im_r_ind({im_bw,r_bw});
	  const qprop_t &prop_inv_in=jprop_inv_in[ip_in][ijack];
	  const qprop_t &prop_inv_ou=jprop_inv_ou[ip_ou][ijack];
	  
	  const qprop_t &vert=jverts[im_r_im_r_iG][ijack];
	  
	  qprop_t amp_vert=prop_inv_ou*vert*quaGamma[5]*prop_inv_in.adjoint()*quaGamma[5];
	  auto c=(amp_vert*quaGamma[iG].adjoint()).trace().real()/(12.0*iG_of_Zbil[iZbil].size());
	  
	  if(0)
	    cout<<"Amputating bil "<<iZbil<<", ip_in: "<<ip_in<<", ip_ou: "<<ip_ou<<", ijack: "<<ijack<<", prop_inv_in: "<<prop_inv_in(0,0)<<", vert: "<<vert(0,0)<<", prop_inv_ou: "<<prop_inv_ou(0,0)<<", c: "<<c<<endl;
	  out[ijack]+=c;
	}
    }
  
  return pr;
}
