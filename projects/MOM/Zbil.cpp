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

void build_jackknifed_vert_Gamma(jprop_t &jvert,const prop_t &prop1,size_t iG,const prop_t &prop2,double w,size_t ijack)
{jvert[ijack]+=w*prop1*Gamma[iG]*Gamma[5]*prop2.adjoint()*Gamma[5];}

void build_all_mr_gbil_jackknifed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_props_t> &props1,const vector<m_r_mom_conf_props_t> &props2,
					const index_t &im_r_im_r_igam_ind,const index_t &im_r_ijack_ind,bool use_QED,const djvec_t &deltam_cr)
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
      const m_r_mom_conf_props_t &p1=props1[im_r_ijack_ind({im_fw,r_fw,ijack})];
      const m_r_mom_conf_props_t &p2=props2[im_r_ijack_ind({im_bw,r_bw,ijack})];
      
      //create list of operations
      vector<tuple<jprop_t*,const prop_t*,const prop_t*,double>> list={{&jbil.LO[im_r_im_r_igam],&p1.LO,&p2.LO,1.0}};
      if(use_QED)
	for(auto &o : vector<tuple<jprop_t*,const prop_t*,const prop_t*,double>>({
	      {&jbil.EM[im_r_im_r_igam],&p1.F,&p2.F,1.0},
	      {&jbil.EM[im_r_im_r_igam],&p1.FF,&p2.LO,1.0},
	      {&jbil.EM[im_r_im_r_igam],&p1.LO,&p2.FF,1.0},
	      {&jbil.EM[im_r_im_r_igam],&p1.T,&p2.LO,1.0},
	      {&jbil.EM[im_r_im_r_igam],&p1.LO,&p2.T,1.0},
	      {&jbil.EM[im_r_im_r_igam],&p1.P,&p2.LO,-deltam_cr[im_fw][ijack]},
	      {&jbil.EM[im_r_im_r_igam],&p1.LO,&p2.P,-deltam_cr[im_bw][ijack]},
	      {&jbil.S[im_r_im_r_igam],&p1.S,&p2.LO,1.0},
	      {&jbil.S[im_r_im_r_igam],&p1.LO,&p2.S,1.0}}))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackknifed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),get<3>(o),ijack);
    }
}

djvec_t compute_proj_bil(const vjprop_t &jprop_inv1,const vector<jprop_t> &jverts,const vjprop_t &jprop_inv2,const index_t &im_r_ind)
{
  const size_t nm=im_r_ind.max(0),nr=im_r_ind.max(1);
  const index_t im_r_im_r_iG_ind({{"im",nm},{"r",nr},{"im",nm},{"r",nr},{"igamma",nGamma}});
  const index_t im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"Zbil",nZbil}});
  const index_t ind({{"rest",im_r_im_r_iZbil_ind.max()},{"ijack",njacks+1}});
  
  //combine the 16 bilinears to form the 6 vertex
  const vector<vector<size_t>> iG_of_Zbil={{0},{1,2,3,4},{5},{6,7,8,9},{10,11,12,13,14,15}};
  djvec_t pr(im_r_im_r_iZbil_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_ibil and ijack
      vector<size_t> i_comp=ind(i);
      const size_t im_r_im_r_ibil=i_comp[0],ijack=i_comp[1];
      const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_iZbil_ind(im_r_im_r_ibil);
      
      //get im and r for fw and back, and ibil
      const size_t im_fw=im_r_im_r_ibil_comp[0],r_fw=im_r_im_r_ibil_comp[1];
      const size_t im_bw=im_r_im_r_ibil_comp[2],r_bw=im_r_im_r_ibil_comp[3];
      const size_t ibil=im_r_im_r_ibil_comp[4];
      
      //loop on all gammas
      djack_t &out=pr[im_r_im_r_ibil];
      for(auto & iG : iG_of_Zbil[ibil])
	{
	  vector<size_t> im_r_im_r_iG_comp=im_r_im_r_ibil_comp;
	  im_r_im_r_iG_comp[4]=iG;
	  const size_t im_r_im_r_iG=im_r_im_r_iG_ind(im_r_im_r_iG_comp);
	  
	  const size_t ip1=im_r_ind({im_fw,r_fw});
	  const size_t ip2=im_r_ind({im_bw,r_bw});
	  const prop_t &prop_inv1=jprop_inv1[ip1][ijack];
	  const prop_t &prop_inv2=jprop_inv2[ip2][ijack];
	  
	  const prop_t &vert=jverts[im_r_im_r_iG][ijack];
	  
	  prop_t amp_vert=prop_inv1*vert*Gamma[5]*prop_inv2.adjoint()*Gamma[5];
	  out[ijack]+=(amp_vert*Gamma[iG].adjoint()).trace().real()/(12.0*iG_of_Zbil[ibil].size());
	}
    }
  
  return pr;
}
