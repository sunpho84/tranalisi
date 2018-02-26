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
      
      auto op=lepGamma[iGl]*(lepGamma[0]
 			     // *0.0
// #warning avoiding identity in mesloop
			     +
			     // 0.0*
			     // sign*
			     sign*lepGamma[5]);
      auto pF=lepGamma[5]*pl.F.adjoint()*lepGamma[5];
      auto pr=(lepGamma[ipGl]*(lepGamma[0]
			       // *0.0
// #warning avoiding identity in mesloop
			       +
			       // 0.0*
			       // psign*
			       psign*lepGamma[5])).adjoint()
// #warning suppressing fact 1/2 in mesloop
	/2.0
	;
      mesloop.LO[i]=(op*pr).toDense().trace()/4.0; //for test: this must be 1 if iGl==ipGl
      mesloop.F[i]=(op*pF*pr).trace()*(double)V/4.0; //normalization for the single gamma
    }
  
  return mesloop;
}

void build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props1,const vector<m_r_mom_conf_qprops_t> &props2,
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
  index_t im_r_im_r_iop_ilistpGl_ind({{"m_fw",nm},{"r_fw",nr},{"m_bw",nm},{"r_bw",nr},{"iop",nZop},{"listpGl",listpGl.size()}});
  
#pragma omp parallel for
   for(size_t i=0;i<im_r_im_r_iop_ilistpGl_iclust_ind.max();i++)
     {
       //decript bilinar/jack
       vector<size_t> comps=im_r_im_r_iop_ilistpGl_iclust_ind(i);
       
       //decript props to combine and which Gamma to create
       const size_t im_fw=comps[0],r_fw=comps[1];
       const size_t im_bw=comps[2],r_bw=comps[3];
       const size_t iop=comps[4],ilistpGl=comps[5];
       const size_t iclust=comps[6];
       
       //proxy for vector and props
       const m_r_mom_conf_qprops_t &p1=props1[im_r_iclust_ind({im_fw,r_fw,iclust})];
       const m_r_mom_conf_qprops_t &p2=props2[im_r_iclust_ind({im_bw,r_bw,iclust})];
       
       //create list of operations
       vector<tuple<vector<jqprop_t>*,const vector<dcompl_t>*,const qprop_t*,const qprop_t*>> list={
	 {&j.LO,&mesloop.LO,&p1.LO,&p2.LO}, //LO
	 
	 
	 //{&j.PH,&mesloop.F,&p1.F,&p2.LO}, //nasty1
	  {&j.PH,&mesloop.F,&p1.LO,&p2.F}, //nasty2
	 
	 // {&j.PH,&mesloop.LO,&p1.FF,&p2.LO}, //self1
	 // {&j.PH,&mesloop.LO,&p1.LO,&p2.FF}, //self2
	 // //
	 // {&j.PH,&mesloop.LO,&p1.T,&p2.LO}, //tad1
	 // {&j.PH,&mesloop.LO,&p1.LO,&p2.T}, //tad2
	 //
	 //{&j.PH,&mesloop.LO,&p1.F,&p2.F}, //exchange
	 //
	 //
	 {&j.CT1,&mesloop.LO,&p1.P,&p2.LO}, //counterterm1
	 {&j.CT2,&mesloop.LO,&p1.LO,&p2.P}, //counterterm2
	 
	 
	 {&j.S,&mesloop.LO,&p1.S,&p2.LO}, //mass1
	 {&j.S,&mesloop.LO,&p1.LO,&p2.S}  //mass2
	 };
       
       const Zop_t &zop=zops[iop];
       
       //create the vertex
       const size_t im_r_im_r_iop_ilistpGl=im_r_im_r_iop_ilistpGl_ind({im_fw,r_fw,im_bw,r_bw,iop,ilistpGl});
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
	       const qprop_t &prop1=*get<2>(o);
	       const qprop_t &prop2=*get<3>(o);
	       
	       const qprop_t c=prop1*quaGamma[Gq]*(quaGamma[0]
 						   // *0.0
// #warning avoiding inserting gamma5 op on quark side
						   +
						   // 0.0*
						   // sign*
						   sign*quaGamma[5])*quaGamma[5]*prop2.adjoint()*quaGamma[5]*mesloop;
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
	       	 ", prop1: "<<prop1(0,0)<<
	       	 ", prop2: "<<prop2(0,0)<<
	       	 ", res: "<<c(0,0)<<
	       	 endl;
	     }
	 }
     }
}

djvec_t compute_proj_meslep(const vjqprop_t &jprop_inv1,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv2,const index_t &im_r_ind)
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
      
      //get im and r for fw and back, and iop and iproj
      const size_t im_fw=im_r_im_r_iop_iproj_comps[0],r_fw=im_r_im_r_iop_iproj_comps[1];
      const size_t im_bw=im_r_im_r_iop_iproj_comps[2],r_bw=im_r_im_r_iop_iproj_comps[3];
      const size_t iop=im_r_im_r_iop_iproj_comps[4],iproj=im_r_im_r_iop_iproj_comps[5];
      
      const double norm=12.0*zops[iproj].norm
	//#warning suppressing factor 1/2 in projection normalization
	*2.0 //2 comes form 1-g5 normalization
	;
      const int pQg5_sign=zops[iproj].Qg5_sign; //same sign
      const size_t ip1=im_r_ind({im_fw,r_fw});
      const size_t ip2=im_r_ind({im_bw,r_bw});
      const qprop_t &prop_inv1=jprop_inv1[ip1][ijack];
      const qprop_t &prop_inv2=jprop_inv2[ip2][ijack];
      
      //loop on all operators
      double &out=pr[im_r_im_r_iop_iproj][ijack];
      out=0.0;
      
      for(auto &pcontr : zops[iproj].contr) // projecting for operator is stored in contr
	{
	  const size_t ilistpGl=pcontr.ilistGl;
	  const size_t im_r_im_r_iop_ilistpGl=im_r_im_r_iop_ilistpGl_ind({im_fw,r_fw,im_bw,r_bw,iop,ilistpGl});
	  const qprop_t &vert=jverts[im_r_im_r_iop_ilistpGl][ijack];
	  
	  const qprop_t amp_vert=prop_inv1*vert*quaGamma[5]*prop_inv2.adjoint()*quaGamma[5];
	  
	  //projecting on quark side
	  const size_t Gq=pcontr.Gq;
	  auto projector=(quaGamma[Gq]*(quaGamma[0]
 					// *0.0
// #warning avoiding gamma5 projecting on quark side
					+
					//0.0*
					// pQg5_sign*
					pQg5_sign*quaGamma[5])).adjoint();
	  auto contr=(amp_vert*projector).trace().real()/norm;
	  out+=contr;
	  
	  if(0)
	  cout<<" amputating, i: "<<i<<", iop: "<<iop<<", ilistpGl: "<<ilistpGl<<", iproj: "<<iproj<<", ijack: "<<ijack<<
	    ", im_r_im_r_iop_ilistpGl(jverts): "<<im_r_im_r_iop_ilistpGl<<", im_r_im_r_iop_iproj(out): "<<im_r_im_r_iop_iproj<<
	    ", prop_inv1: "<<prop_inv1(0,0)<<", vert: "<<vert(0,0)<<", prop_inv2: "<<prop_inv2(0,0)<<", contr: "<<contr<<endl;
	  }
      
      // cout<<"amputated, im_fw: "<<im_fw<<", r_fw: "<<r_fw<<", im_bw: "<<im_bw<<", r_bw: "<<r_bw<<", iop: "<<iop<<", iproj: "<<iproj<<", ijack: "<<ijack<<", out: "<<out<<endl;
    }
  
  return pr;
}
