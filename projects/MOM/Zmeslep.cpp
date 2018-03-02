#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_MESLEP
 #include <Zmeslep.hpp>

#include <geometry.hpp>
#include <ingredients.hpp>
#include <oper.hpp>

#include <prop.hpp>

namespace meslep
{
  index_t get_ilistGl_ilistpGl_iclust_ind()
  {
    return index_t({{"listGl",listGl.size()},{"pGl",listpGl.size()},{"iclust",njacks}});
  }
}

meslep::mesloop_t build_mesloop(const vector<mom_conf_lprops_t> &props_lep)
{
  using namespace meslep;
  
  const index_t ilistGl_ilistpGl_iclust_ind=get_ilistGl_ilistpGl_iclust_ind();
  
  //! projected lepton propagator with insertion
  meslep::mesloop_t mesloop(ilistGl_ilistpGl_iclust_ind.max());
  
  // NB: the lepton loop is fully amputated
#pragma omp parallel for
  for(size_t i=0;i<ilistGl_ilistpGl_iclust_ind.max();i++)
    {
      vector<size_t> comps=ilistGl_ilistpGl_iclust_ind(i);
      const size_t ilistGl=comps[0],ilistpGl=comps[1],iclust=comps[2];
      const mom_conf_lprops_t &pl=props_lep[iclust];
      const size_t iGl=listGl[ilistGl];
      const size_t ipGl=listpGl[ilistpGl];
      const int sign=Lg5_sign[ilistGl];
      const int psign=Lg5_sign[ilistpGl];
      
      auto op=lepGamma[iGl]*(lepGamma[0]+sign*lepGamma[5]);
      lprop_t p0=lepGamma[5]*pl.LO.adjoint()*lepGamma[5];
      lprop_t pF=lepGamma[5]*pl.F.adjoint()*lepGamma[5];
      lprop_t pFamp=pF*p0.inverse();
      auto pr=(lepGamma[ipGl]*(lepGamma[0]+psign*lepGamma[5])).adjoint()/2.0;
      mesloop.LO[i]=(op*pr).toDense().trace()/4.0;   //for test: this must be 1 if iGl==ipGl
      mesloop.F[i]=(op*pFamp*pr).trace()*(double)V/4.0; //normalization for the single gamma
    }
  
  return mesloop;
}

void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou,
					    const vector<mom_conf_lprops_t> &props_lep,const index_t &im_r_ind)
{
  using namespace meslep;
  
  const size_t nm=im_r_ind.max(0);
  const size_t nr=im_r_ind.max(1);
  
  const mesloop_t mesloop=build_mesloop(props_lep);
  
  const index_t ilistGl_ilistpGl_iclust_ind=get_ilistGl_ilistpGl_iclust_ind();
  
  //! help finding the meslep/clust combo
  index_t im_r_iclust_ind=im_r_ind*index_t({{"clust",njacks}});
  index_t im_r_im_r_iop_ilistpGl_iclust_ind=im_r_ind*im_r_ind*index_t({{"iop",nZop},{"listpGl",listpGl.size()},{"clust",njacks}});
  index_t im_r_im_r_iop_ilistpGl_ind({{"m_in",nm},{"r_in",nr},{"m_ou",nm},{"r_ou",nr},{"iop",nZop},{"listpGl",listpGl.size()}});
  
#pragma omp parallel for
   for(size_t i=0;i<im_r_im_r_iop_ilistpGl_iclust_ind.max();i++)
     {
       //decript bilinar/jack
       vector<size_t> comps=im_r_im_r_iop_ilistpGl_iclust_ind(i);
       
       //decript props to combine and which Gamma to create
       const size_t im_in=comps[0],r_in=comps[1];
       const size_t im_ou=comps[2],r_ou=comps[3];
       const size_t iop=comps[4],ilistpGl=comps[5];
       const size_t iclust=comps[6];
       
       //proxy for vector and props
       const m_r_mom_conf_qprops_t &p_in=props_in[im_r_iclust_ind({im_in,r_in,iclust})];
       const m_r_mom_conf_qprops_t &p_ou=props_ou[im_r_iclust_ind({im_ou,r_ou,iclust})];
       
       //create list of operations
       vector<tuple<vector<jqprop_t>*,const vector<dcompl_t>*,const qprop_t*,const qprop_t*>> list={
	 {&j.LO,&mesloop.LO,&p_in.LO,&p_ou.LO}, //LO
	 //
	 //
	 {&j.PH,&mesloop.F,&p_in.F,&p_ou.LO}, //nasty_in
	 {&j.PH,&mesloop.F,&p_in.LO,&p_ou.F}, //nasty_ou
	 //
	 {&j.PH,&mesloop.LO,&p_in.FF,&p_ou.LO}, //self_in
	 {&j.PH,&mesloop.LO,&p_in.LO,&p_ou.FF}, //self_ou
	 //
	 {&j.PH,&mesloop.LO,&p_in.T,&p_ou.LO}, //tad_in
	 {&j.PH,&mesloop.LO,&p_in.LO,&p_ou.T}, //tad_ou
	 //
	 {&j.PH,&mesloop.LO,&p_in.F,&p_ou.F}, //exchange
	 //
	 //
	 {&j.CT_IN,&mesloop.LO,&p_in.P,&p_ou.LO}, //counterterm_in
	 {&j.CT_OU,&mesloop.LO,&p_in.LO,&p_ou.P}, //counterterm_ou
	 //
	 //
	 {&j.S,&mesloop.LO,&p_in.S,&p_ou.LO}, //mass_in
	 {&j.S,&mesloop.LO,&p_in.LO,&p_ou.S}  //mass_ou
	 };
       
       const Zop_t &zop=zops[iop];
       
       //create the vertex
       const size_t im_r_im_r_iop_ilistpGl=im_r_im_r_iop_ilistpGl_ind({im_in,r_in,im_ou,r_ou,iop,ilistpGl});
       for(auto &o : list)
	 {
	   const int    sign=zop.Qg5_sign;
	   jqprop_t &jvert=(*get<0>(o))[im_r_im_r_iop_ilistpGl];
	   // cout<<" "<<endl;
	   for(auto &contr : zop.contr)
	     {
	       const size_t ilistGl=contr.ilistGl;
	       const size_t Gq=contr.Gq;
	       const size_t imesloop=ilistGl_ilistpGl_iclust_ind({ilistGl,ilistpGl,iclust});
	       const dcompl_t &mesloop=(*get<1>(o))[imesloop];
	       const qprop_t &prop_in=*get<2>(o);
	       const qprop_t &prop_ou=*get<3>(o);
	       
	       const qprop_t c=prop_ou*quaGamma[Gq]*(quaGamma[0]+sign*quaGamma[5])*quaGamma[5]*prop_in.adjoint()*quaGamma[5]*mesloop;
	       jvert[iclust]+=c;
	       
	       if(0)
	       cout
	       	 <<"iop: "<<iop<<
	       	 ", ilistGl: "<<ilistGl<<"("<<listGl[ilistGl]<<")"<<
	       	 ", Gq: "<<Gq<<
	       	 ", sign: "<<sign<<
	       	 ", ilistpGl: "<<ilistpGl<<"("<<listpGl[ilistpGl]<<")"<<
	       	 ", iclust: "<<iclust<<
	       	 ", mesloop: "<<mesloop<<
	       	 ", prop_in: "<<prop_in(0,0)<<
	       	 ", prop_ou: "<<prop_ou(0,0)<<
	       	 ", res: "<<c(0,0)<<
	       	 endl;
	     }
	 }
     }
}

djvec_t compute_proj_meslep(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou,const index_t &im_r_ind)
{
  using namespace meslep;
  
  const index_t im_r_im_r_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"iop",nZop},{"iproj",nZop}});
  const index_t im_r_im_r_iop_ilistpGl_ind=im_r_ind*im_r_ind*index_t({{"iop",nZop},{"listpGL",listpGl.size()}});
  const index_t ind({{"rest",im_r_im_r_iop_iproj_ind.max()},{"ijack",njacks+1}});
  
  djvec_t pr(im_r_im_r_iop_iproj_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<ind.max();i++)
    {
      //split im_r_im_r_iop_iproj and ijack
      vector<size_t> comps=ind(i);
      const size_t im_r_im_r_iop_iproj=comps[0],ijack=comps[1];
      const vector<size_t> im_r_im_r_iop_iproj_comps=im_r_im_r_iop_iproj_ind(im_r_im_r_iop_iproj);
      
      //get im and r for in and back, and iop and iproj
      const size_t im_in=im_r_im_r_iop_iproj_comps[0],r_in=im_r_im_r_iop_iproj_comps[1];
      const size_t im_ou=im_r_im_r_iop_iproj_comps[2],r_ou=im_r_im_r_iop_iproj_comps[3];
      const size_t iop=im_r_im_r_iop_iproj_comps[4],iproj=im_r_im_r_iop_iproj_comps[5];
      
      const double norm=12.0*zops[iproj].norm*2.0;
      const int pQg5_sign=zops[iproj].Qg5_sign; //same sign
      const size_t ip_in=im_r_ind({im_in,r_in});
      const size_t ip_ou=im_r_ind({im_ou,r_ou});
      const qprop_t &prop_inv_in=jprop_inv_in[ip_in][ijack];
      const qprop_t &prop_inv_ou=jprop_inv_ou[ip_ou][ijack];
      
      //loop on all operators
      double &out=pr[im_r_im_r_iop_iproj][ijack];
      out=0.0;
      
      for(auto &pcontr : zops[iproj].contr) // projecting for operator is stored in contr
	{
	  const size_t ilistpGl=pcontr.ilistGl;
	  const size_t im_r_im_r_iop_ilistpGl=im_r_im_r_iop_ilistpGl_ind({im_in,r_in,im_ou,r_ou,iop,ilistpGl});
	  const qprop_t &vert=jverts[im_r_im_r_iop_ilistpGl][ijack];
	  
	  const qprop_t amp_vert=prop_inv_ou*vert*quaGamma[5]*prop_inv_in.adjoint()*quaGamma[5];
	  
	  //projecting on quark side
	  const size_t Gq=pcontr.Gq;
	  auto projector=(quaGamma[Gq]*(quaGamma[0]+pQg5_sign*quaGamma[5])).adjoint();
	  auto contr=(amp_vert*projector).trace().real()/norm;
	  out+=contr;
	  
	  if(0)
	  cout<<" amputating, i: "<<i<<", iop: "<<iop<<", ilistpGl: "<<ilistpGl<<", iproj: "<<iproj<<", ijack: "<<ijack<<
	    ", im_r_im_r_iop_ilistpGl(jverts): "<<im_r_im_r_iop_ilistpGl<<", im_r_im_r_iop_iproj(out): "<<im_r_im_r_iop_iproj<<
	    ", prop_inv_in: "<<prop_inv_in(0,0)<<", vert: "<<vert(0,0)<<", prop_inv_ou: "<<prop_inv_ou(0,0)<<", contr: "<<contr<<endl;
	  }
      
      // cout<<"amputated, im_in: "<<im_in<<", r_in: "<<r_in<<", im_ou: "<<im_ou<<", r_ou: "<<r_ou<<", iop: "<<iop<<", iproj: "<<iproj<<", ijack: "<<ijack<<", out: "<<out<<endl;
    }
  
  return pr;
}
