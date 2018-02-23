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

void build_jackkniffed_vert_Gamma(jqprop_t &jvert,const qprop_t &prop1,size_t iG,const qprop_t &prop2,size_t iclust)
{
  auto c=prop1*quaGamma[iG]*quaGamma[5]*prop2.adjoint()*quaGamma[5];
  jvert[iclust]+=c;
  
  cout
    <<"iG: "<<iG<<
    ", iclust: "<<iclust<<
    ", prop1: "<<prop1(0,0)<<
    ", prop2: "<<prop2(0,0)<<
    ", res: "<<c(0,0)<<
    endl;

}

void build_all_mr_gbil_jackkniffed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
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
      const m_r_mom_conf_qprops_t &p1=props1[im_r_iclust_ind({im_fw,r_fw,iclust})];
      const m_r_mom_conf_qprops_t &p2=props2[im_r_iclust_ind({im_bw,r_bw,iclust})];
      
      //create list of operations
      vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>> list


	// ={{&jbil.LO[im_r_im_r_igam],&p1.LO,&p2.LO}}


	;
      if(use_QED)
	for(auto &o : vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>>({
	      {&jbil.PH[im_r_im_r_igam],&p1.F,&p2.F},
	      // {&jbil.PH[im_r_im_r_igam],&p1.FF,&p2.LO},
	      // {&jbil.PH[im_r_im_r_igam],&p1.LO,&p2.FF},
	      // {&jbil.PH[im_r_im_r_igam],&p1.T,&p2.LO},
	      // {&jbil.PH[im_r_im_r_igam],&p1.LO,&p2.T},
	      // {&jbil.CT1[im_r_im_r_igam],&p1.P,&p2.LO},
	      // {&jbil.CT2[im_r_im_r_igam],&p1.LO,&p2.P},
	      // {&jbil.S[im_r_im_r_igam],&p1.S,&p2.LO},
	      // {&jbil.S[im_r_im_r_igam],&p1.LO,&p2.S}
		}))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackkniffed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),iclust);
    }
}

djvec_t compute_proj_bil(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind)
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
	  
	  const size_t ip1=im_r_ind({im_fw,r_fw});
	  const size_t ip2=im_r_ind({im_bw,r_bw});
	  const qprop_t &prop_inv1=jprop_inv1[ip1][ijack];
	  const qprop_t &prop_inv2=jprop_inv2[ip2][ijack];
	  
	  const qprop_t &vert=jverts[im_r_im_r_iG][ijack];
	  
	  qprop_t amp_vert=prop_inv1*vert*quaGamma[5]*prop_inv2.adjoint()*quaGamma[5];
	  auto c=(amp_vert*quaGamma[iG].adjoint()).trace().real()/(12.0*iG_of_Zbil[iZbil].size());
	  
	  cout<<"Amputating bil "<<iZbil<<", ip1: "<<ip1<<", ip2: "<<ip2<<", prop_inv1: "<<prop_inv1(0,0)<<", prop_inv2: "<<prop_inv2(0,0)<<", c: "<<c<<endl;
	  out[ijack]+=c;
	}
    }
  
  return pr;
}
