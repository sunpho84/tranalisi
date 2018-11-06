#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_PR_MESLEP
 #include <MOM2/pr_meslep.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

using namespace meslep;

namespace pr_meslep
{
  void set_ins()
  {
    switch(pars::use_QED)
      {
      case 0:
      	ins_list={LO};
	break;
      case 1:
	ins_list={LO , QED , NA_IN , NA_OU , CR_OU , CR_IN , TM_OU , TM_IN , PH_OU , PH_IN , EX };
	break;
      case 2:
	ins_list={LO , QED , NA_IN , NA_OU , QED_OU , QED_IN , EX };
	break;
      }
    
    iins_of_ins.resize(ins_tag.size());
    for(size_t iins=0;iins<ins_list.size();iins++)
      iins_of_ins[ins_list[iins]]=iins;
    nins=ins_list.size();
    cout<<"Pr_meslep, nins: "<<nins<<endl;
  }
 }

vector<dcompl_t> perens_t::build_mesloop(const vector<lprop_t> &props_lep) const
{
  //! projected lepton propagator with insertion
  vector<dcompl_t> mesloop(ilistGl_ilistpGl_lins_iclust_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<ilistGl_ilistpGl_iclust_ind.max();i++)
    {
      vector<size_t> comps=ilistGl_ilistpGl_iclust_ind(i);
      const size_t ilistGl=comps[0],ilistpGl=comps[1],iclust=comps[2];
      const size_t iGl=listGl[ilistGl];
      const size_t ipGl=listpGl[ilistpGl];
      const int sign=Lg5_sign[ilistGl];
      const int psign=Lg5_sign[ilistpGl];
      
      const auto op=lepGamma[iGl]*(lepGamma[0]+sign*lepGamma[5]);
      const auto p0=lepGamma[5]*props_lep[ilins_ijack_ind({lprop::LO,iclust})].adjoint()*lepGamma[5];
      const auto pF=lepGamma[5]*props_lep[ilins_ijack_ind({lprop::F,iclust})].adjoint()*lepGamma[5];
      const auto pFamp=pF*p0.inverse();
      const auto pr=(lepGamma[ipGl]*(lepGamma[0]+psign*lepGamma[5])).adjoint()/2.0;
      auto ml=[&](lprop::ins ins)->dcompl_t&{return mesloop[ilistGl_ilistpGl_lins_iclust_ind({{ilistGl,ilistpGl,ins,iclust}})];};
      ml(lprop::LO)=(op*pr).toDense().trace()/4.0;   //for test: this must be 1 if iGl==ipGl
      ml(lprop::F)=(op*pFamp*pr).trace()/4.0;        //normalization for the single gamma
      //cout<<"p0 "<<p0(0,0)<<"\npF "<<pF(0,0)<<"\npFamp "<<pFamp(0,0)<<"\npr "<<((lprop_t)pr)(0,0)<<"\nml "<<ml(lprop::LO)<<"\nmlF "<<ml(lprop::F)<<endl;
    }
  
  return mesloop;
}

void perens_t::build_all_mr_gmeslep_jackkniffed_verts(vector<jqprop_t> &j,const vector<qprop_t> &props_in,const vector<qprop_t> &props_ou,
						      const vector<lprop_t> &props_lep) const
{
  const vector<dcompl_t> mesloop=build_mesloop(props_lep);
  
  //! help finding the meslep/clust combo
  index_t im_r_im_r_iop_ilistpGl_iclust_ind=im_r_ind*im_r_ind*index_t({{"iop",nZop},{"listpGl",listpGl.size()},{"clust",njacks}});
  index_t im_r_im_r_iop_ilistpGl_ind({{"m_in",nm},{"r_in",nr},{"m_ou",nm},{"r_ou",nr},{"iop",nZop},{"listpGl",listpGl.size()}});
  
  //! list of all combination of transformations to be applied
  vector<tuple<size_t,lprop::ins,size_t,size_t>> map;
  
#define ADD_COMBO(P,L,I,O)					\
  map.push_back({pr_meslep::iins_of_ins[pr_meslep::P],lprop::L,qprop::iins_of_ins[qprop::I],qprop::iins_of_ins[qprop::O]})
  
  ADD_COMBO(    LO,    LO, LO, LO); //leading order
  
  switch(pars::use_QED)
    {
    case 0:
      break;
    case 1:
      ADD_COMBO(EX,    LO, F,F);   //exchange
      //
      ADD_COMBO(NA_IN, F,  F,LO);  //nasty_in
      ADD_COMBO(NA_OU, F,  LO,F);  //nasty_ou
      //
      ADD_COMBO(PH_IN, LO, FF,LO); //self_in
      ADD_COMBO(PH_OU, LO, LO,FF); //self_ou
      //
      ADD_COMBO(PH_IN, LO, T,LO);  //tad_in
      ADD_COMBO(PH_OU, LO, LO,T);  //tad_ou
      //
      ADD_COMBO(CR_IN, LO, P,LO);  //critical counterterm_in
      ADD_COMBO(CR_OU, LO, LO,P);  //critical counterterm_ou
      //
      ADD_COMBO(TM_IN, LO, S,LO);  //twisted counterterm_in
      ADD_COMBO(TM_OU, LO, LO,S);  //twisted counterterm_ou
      break;
    case 2:
      ADD_COMBO(EX,    LO, F,F);   //exchange
      //
      ADD_COMBO(NA_IN, F,  F,LO);  //nasty_in
      ADD_COMBO(NA_OU, F,  LO,F);  //nasty_ou
      //
      ADD_COMBO(QED_IN, LO, QED,LO); //all QED_in
      ADD_COMBO(QED_OU, LO, LO,QED); //all QED_ou
      break;
    }
#undef ADD_COMBO
  
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
       
       const Zop_t &zop=zops[iop];
       
       //create the vertex
       for(auto &o : map)
	 {
	   const size_t v_ins=get<0>(o);
	   const lprop::ins ml_ins=get<1>(o);
	   const size_t iq_ins_in=get<2>(o);
	   const size_t iq_ins_ou=get<3>(o);
	   
	   const int    sign=zop.Qg5_sign;
	   const int    norm=zop.norm;
	   
	   qprop_t &vert=j[im_r_im_r_iop_ilistpGl_meslepins_ind({im_in,r_in,im_ou,r_ou,iop,ilistpGl,v_ins})][iclust];
	   for(auto &contr : zop.contr)
	     {
	       const size_t ilistGl=contr.ilistGl;
	       const size_t Gq=contr.Gq;
	       const size_t imesloop=ilistGl_ilistpGl_lins_iclust_ind({{ilistGl,ilistpGl,ml_ins,iclust}});
	       const dcompl_t &ml=mesloop[imesloop];
	       const qprop_t &prop_in=props_in[im_r_iqins_ijack_ind({im_in,r_in,iq_ins_in,iclust})]; ;
	       const qprop_t &prop_ou=props_ou[im_r_iqins_ijack_ind({im_ou,r_ou,iq_ins_ou,iclust})]; ;
	       
	       const qprop_t c=prop_ou*quaGamma[Gq]*(quaGamma[0]+sign*quaGamma[5])*quaGamma[5]*prop_in.adjoint()*quaGamma[5]*ml;
	       vert+=c*norm;
	       
	       if(0)
		 cout
		   <<"iop: "<<iop<<
		   ", ilistGl: "<<ilistGl<<"("<<listGl[ilistGl]<<")"<<
		   ", Gq: "<<Gq<<
		   ", sign: "<<sign<<
		   ", ilistpGl: "<<ilistpGl<<"("<<listpGl[ilistpGl]<<")"<<
		   ", iclust: "<<iclust<<
		   ", mesloop: "<<ml<<
		   ", prop_in: "<<prop_in(0,0)<<
		   ", prop_ou: "<<prop_ou(0,0)<<
		   ", res: "<<c(0,0)<<
		   endl;
	     }
	   
	   //cout<<" "<<vert<<endl;
	 }
     }
}

void perens_t::compute_proj_meslep(const vector<jqprop_t> &jprop_inv_in,const vector<jqprop_t> &jverts,const vector<jqprop_t> &jprop_inv_ou,const size_t imeslepmom)
{
  const index_t ind({{"rest",im_r_im_r_iop_iproj_ind.max()},{"ijack",njacks+1}});
  
  vector<tuple<size_t,size_t,size_t,size_t>> map;
#define ADD_COMBO(A,I,V,O)						\
  map.push_back({pr_meslep::iins_of_ins[pr_meslep::A],jqprop::iins_of_ins[jqprop::I],pr_meslep::iins_of_ins[pr_meslep::V],jqprop::iins_of_ins[jqprop::O]})
  ADD_COMBO(LO,  LO, LO, LO);
  switch(pars::use_QED)
    {
    case 0:
      break;
    case 1:
      ADD_COMBO(EX    , LO , EX    , LO);
      ADD_COMBO(NA_IN , LO , NA_IN , LO);
      ADD_COMBO(NA_OU , LO , NA_OU , LO);
      ADD_COMBO(CR_IN , LO , CR_IN , LO);
      ADD_COMBO(CR_OU , LO , CR_OU , LO);
      ADD_COMBO(TM_IN , LO , TM_IN , LO);
      ADD_COMBO(TM_OU , LO , TM_OU , LO);
      ADD_COMBO(PH_IN , LO , PH_IN , LO);
      ADD_COMBO(PH_OU , LO , PH_OU , LO);
      
      ADD_COMBO(CR_IN , CR , LO, LO);
      ADD_COMBO(CR_OU , LO , LO, CR);
      ADD_COMBO(TM_IN , TM , LO, LO);
      ADD_COMBO(TM_OU , LO , LO, TM);
      ADD_COMBO(PH_IN , PH , LO, LO);
      ADD_COMBO(PH_OU , LO , LO, PH);
      break;
    case 2:
      ADD_COMBO(EX     , LO  , EX     , LO);
      ADD_COMBO(NA_IN  , LO  , NA_IN  , LO);
      ADD_COMBO(NA_OU  , LO  , NA_OU  , LO);
      ADD_COMBO(QED_IN , LO  , QED_IN , LO);
      ADD_COMBO(QED_OU , LO  , QED_OU , LO);
      
      ADD_COMBO(QED_IN , QED , LO     , LO);
      ADD_COMBO(QED_OU , LO  , LO     , QED);
      break;
    }
#undef ADD_COMBO

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
      
      const double pnorm=12.0*zops[iproj].pnorm*2.0;
      const int pQg5_sign=zops[iproj].Qg5_sign; //same sign
      
      for(auto t : map)
	{
	  const size_t pr_meslepins=get<0>(t);
	  const size_t ijqins_in=get<1>(t);
	  const size_t meslepins=get<2>(t);
	  const size_t ijqins_ou=get<3>(t);
	  
	  const size_t ip_in=im_r_ijqins_ind({im_in,r_in,ijqins_in});
	  const size_t ip_ou=im_r_ijqins_ind({im_ou,r_ou,ijqins_ou});
	  
	  const qprop_t &prop_inv_in=jprop_inv_in[ip_in][ijack];
	  const qprop_t &prop_inv_ou=jprop_inv_ou[ip_ou][ijack];
	  
	  //loop on all operators
	  const size_t im_r_im_r_meslepins_iop_iproj_imeslepmom=
	    im_r_im_r_meslepins_iop_iproj_imeslepmom_ind({im_in,r_in,im_ou,r_ou,pr_meslepins,iop,iproj,imeslepmom});
	  double &out=pr_meslep[im_r_im_r_meslepins_iop_iproj_imeslepmom][ijack];
	  
	  for(auto &pcontr : zops[iproj].contr) // projecting for operator is stored in contr
	    {
	      const size_t ilistpGl=pcontr.ilistGl;
	      const qprop_t &vert=jverts[im_r_im_r_iop_ilistpGl_meslepins_ind({im_in,r_in,im_ou,r_ou,iop,ilistpGl,meslepins})][ijack];
	      
	      const qprop_t amp_vert=prop_inv_ou*vert*quaGamma[5]*prop_inv_in.adjoint()*quaGamma[5];
	      
	      //projecting on quark side
	      const size_t Gq=pcontr.Gq;
	      auto projector=(quaGamma[Gq]*(quaGamma[0]+pQg5_sign*quaGamma[5])).adjoint();
	      auto contr=(amp_vert*projector).trace().real()/pnorm;
	      out+=contr;
	    }
	  
	  //cout<<"amputated, im_in: "<<im_in<<", r_in: "<<r_in<<", im_ou: "<<im_ou<<", r_ou: "<<r_ou<<", iop: "<<iop<<", iproj: "<<iproj<<", ijack: "<<ijack<<", out: "<<out<<endl;
	}
    }
}

void perens_t::mom_compute_meslep()
{
  vector<raw_file_t> qfiles=setup_read_all_qprops_mom(conf_list);
  vector<raw_file_t> lfiles=setup_read_all_lprops_mom(conf_list);
  
  for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
    {
      const size_t imom_in=meslepmoms()[imeslepmom][1];
      const size_t imom_ou=meslepmoms()[imeslepmom][2];
      const bool read_ou=(imom_in!=imom_ou);
      
      vector<jqprop_t> jprops_ou(im_r_ijqins_ind.max());                         //!< jackknived props out
      vector<jqprop_t> jprops_in(im_r_ijqins_ind.max());                         //!< jackknived props in
      vector<jqprop_t> jverts(im_r_im_r_iop_ilistpGl_meslepins_ind.max());       //!< jackknived vertex
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
  	    vector<qprop_t> props_in=read_all_qprops_mom(qfiles,i_in_clust_ihit,mom_in);
  	    vector<qprop_t> props_ou=(read_ou?read_all_qprops_mom(qfiles,i_in_clust_ihit,mom_ou):props_in);
  	    vector<lprop_t> props_lep=read_all_lprops_mom(lfiles,i_in_clust_ihit,momlep);
  	    read_time.stop();
	    
  	    //build all props
  	    build_props_time.start();
  	    build_all_mr_jackkniffed_qprops(jprops_in,props_in);
  	    build_all_mr_jackkniffed_qprops(jprops_ou,props_ou);
  	    build_props_time.stop();
	    
  	    //build all meslep_verts
  	    build_meslep_verts_time.start();
  	    build_all_mr_gmeslep_jackkniffed_verts(jverts,props_in,props_ou,props_lep);
  	    build_meslep_verts_time.stop();
  	  }
      
      //clusterize
      clust_time.start();
      for(auto o : {&jprops_ou,&jprops_in,&jverts})
	for(auto &j : *o)
	  j.clusterize(clust_size);
      clust_time.stop();
      
      vector<jqprop_t> jprops_inv_in=get_inverse_propagators(jprops_in); //!< inverse propagator_in
      vector<jqprop_t> jprops_inv_ou=get_inverse_propagators(jprops_ou); //!< inverse propagator_ou
      
      proj_time.start();
      compute_proj_meslep(jprops_inv_in,jverts,jprops_inv_ou,imeslepmom);
      proj_time.stop();
    }
}

void perens_t::subtract_Oa2_pr_meslep()
{
  CRASH("");
}

void perens_t::evolve_pr_meslep(perens_t &out) const
{
  CRASH("");
}

vector<perens_t::task_t> perens_t::get_pr_meslep_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_pr_meslep;
  for(auto &e : ens)
    in_pr_meslep.push_back(&e->pr_meslep);
  
  vector<task_t> pr_meslep_tasks;
  
  if(pars::compute_meslep) pr_meslep_tasks.push_back({&pr_meslep,in_pr_meslep,im_r_im_r_meslepins_iop_iproj_imeslepmom_ind,"pr_meslep",QCD_task});
  
  return pr_meslep_tasks;
}

void perens_t::average_r_pr_meslep(perens_t &out,const array<double,2> w) const
{
  cout<<"Averaging r for pr_meslep"<<endl;
  
  for(auto &t : out.get_pr_meslep_tasks({this}))
    {
      cout<<" "<<t.tag<<endl;
      
      const djvec_t &pr=*t.in.front();
      djvec_t &pr_rave=*t.out;
      
#pragma omp parallel for
      for(size_t out_i=0;out_i<out.im_r_im_r_meslepins_iop_iproj_imeslepmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_iop_iproj_imeslepmom_comp=out.im_r_im_r_meslepins_iop_iproj_imeslepmom_ind(out_i);
	  vector<size_t> im_r_im_r_iop_iproj_imeslepmom_comp=out_im_r_im_r_iop_iproj_imeslepmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_im_r_iop_iproj_imeslepmom_comp[1]=
		im_r_im_r_iop_iproj_imeslepmom_comp[3]=r;
	      const size_t i=im_r_im_r_meslepins_iop_iproj_imeslepmom_ind(im_r_im_r_iop_iproj_imeslepmom_comp);
	      pr_rave[out_i]+=w[r]*pr[i];
	    }
	  pr_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_pr_meslep(perens_t &out,const vector<vector<size_t>> &equiv_meslepmom_combos) const
{
#pragma omp parallel for
  for(size_t i=0;i<out.im_r_im_r_meslepins_iop_iproj_imeslepmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_iop_iproj_imeslepmom_comp=out.im_r_im_r_meslepins_iop_iproj_imeslepmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_iop_iproj_imeslepmom_comp[7];
      
      for(const auto &t : out.get_pr_meslep_tasks({this}))
	{
	  djack_t &ave=(*t.out)[i];
	  ave=0.0;
	  for(const size_t ieq_mom : equiv_meslepmom_combos[out_imom_combo])
	    {
	      vector<size_t> im_r_im_r_iop_iproj_imeslepmom_comp=out_im_r_im_r_iop_iproj_imeslepmom_comp;
	      im_r_im_r_iop_iproj_imeslepmom_comp[7]=ieq_mom;
	      const size_t ieq=im_r_im_r_meslepins_iop_iproj_imeslepmom_ind(im_r_im_r_iop_iproj_imeslepmom_comp);
	      ave+=(*t.in.front())[ieq];
	    }
	  ave/=equiv_meslepmom_combos[out_imom_combo].size();
	}
    }
}

void perens_t::val_chir_extrap_pr_meslep(perens_t &out) const
{
  const index_t r_r_iop_iproj_imeslepmom_ind({{"r",nr},{"r",nr},{"op",meslep::nZop},{"proj",meslep::nZop},{"meslepmoms",meslepmoms().size()}});
  const index_t r_r_meslepins_iop_iproj_imeslepmom_ind({{"r",nr},{"r",nr},{"meslepins",pr_meslep::nins},{"op",meslep::nZop},{"proj",meslep::nZop},{"meslepmoms",meslepmoms().size()}});
  
  for(auto &t : out.get_pr_meslep_tasks({this}))
#pragma omp parallel for
    for(size_t r_r_iop_iproj_imeslepmom=0;r_r_iop_iproj_imeslepmom<r_r_iop_iproj_imeslepmom_ind.max();r_r_iop_iproj_imeslepmom++)
      {
	const vector<size_t> r_r_iop_iproj_imeslepmom_comps=r_r_iop_iproj_imeslepmom_ind(r_r_iop_iproj_imeslepmom);
	const size_t r_in=r_r_iop_iproj_imeslepmom_comps[0];
	const size_t r_ou=r_r_iop_iproj_imeslepmom_comps[1];
	const size_t iop=r_r_iop_iproj_imeslepmom_comps[2];
	const size_t iproj=r_r_iop_iproj_imeslepmom_comps[3];
	const size_t imeslepmom=r_r_iop_iproj_imeslepmom_comps[4];
	
	const djvec_t &pr=*t.in.front();
	djvec_t &pr_chir=*t.out;
	const string &tag=t.tag;
	
	//check if we need to subtract the pole
	const bool sub_pole=(iop==2 or iop==3 or iproj==2 or iproj==3);
	const size_t x_pow=(sub_pole?1:0);
	
	vector<djvec_t> coeffs(pr_meslep::nins);
	for(size_t meslepins=0;meslepins<pr_meslep::nins;meslepins++)
	  {
	    const vector<size_t> r_r_meslepins_iop_iproj_imeslepmom_comps={r_in,r_ou,meslepins,iop,iproj,imeslepmom};
	    const size_t r_r_meslepins_iop_iproj_imeslepmom=r_r_meslepins_iop_iproj_imeslepmom_ind(r_r_meslepins_iop_iproj_imeslepmom_comps);
	    
	    //open the plot file if needed
	    const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+r_r_meslepins_iop_iproj_imeslepmom_ind.descr(r_r_meslepins_iop_iproj_imeslepmom)+".xmg";
	    grace_file_t *plot=nullptr;
	    if(imeslepmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	    
	    //slice m and fit it
	    djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
	    vector<double> x(nm*(nm+1)/2);
	    int i=0;
	    for(size_t im_ou=0;im_ou<nm;im_ou++)
	      for(size_t im_in=im_ou;im_in<nm;im_in++)
		{
		  const size_t imeson=im_im_ind({im_ou,im_in});
		  
		  //compute mass sum
		  if(pars::chir_extr_method==chir_extr::MQUARK) x[i]=am[im_ou]+am[im_in];
		  else                                          x[i]=sqr(meson_mass[imeson].ave());
		  
		  //compute y and y_plot
		  y_plot[i]=pr[im_r_im_r_meslepins_iop_iproj_imeslepmom_ind({im_ou,r_in,im_in,r_ou,meslepins,iop,iproj,imeslepmom})];
		  
		  //if QED case and pole must be subtracted, take into account variation due to leading pole
		  if(pars::use_QED and meslepins>0)
		    {
		      const djack_t M=meson_mass[imeson],dM=meson_mass_QED[imeson];
		      if(sub_pole)
			{
			  const djack_t b0=coeffs[pr_meslep::LO][2],c0=coeffs[pr_meslep::LO][0];
			  const djack_t varb=2.0*b0*dM*M;
			  const djack_t varc=-2.0*c0*dM/(M*M*M);
			  y_plot[i]-=varb+varc;
			}
		      else
			if(pars::sub_meson_mass_shift_when_no_pole)
			  {
			    const djack_t b0=coeffs[pr_meslep::LO][1];
			    const djack_t varb=2.0*b0*dM*M;
			    y_plot[i]-=varb;
			  }
		    }
		  
		  //fit x*y if pole present
		  y[i]=pow(x[i],x_pow)*y_plot[i];
		  
		  i++;
		}
	    
	    //fit, store and write the result
	    coeffs[meslepins]=poly_fit(x,y,(sub_pole?2:1));
	    
	    if(std::isnan(coeffs[meslepins][0][0])) coeffs[meslepins]=0.0;
	    
	    const size_t iout=out.im_r_im_r_meslepins_iop_iproj_imeslepmom_ind({0,r_in,0,r_ou,meslepins,iop,iproj,imeslepmom});
	    pr_chir[iout]=coeffs[meslepins][x_pow];
	    if(plot!=nullptr)
	      {
		auto xminmax=minmax_element(x.begin(),x.end());
		double xmin=*xminmax.first*0.9;
		double xmax=*xminmax.second*1.1;
		write_fit_plot(*plot,xmin,xmax,[&coeffs,sub_pole,meslepins](double x)->djack_t{return poly_eval<djvec_t>(coeffs[meslepins],x)/(sub_pole?x:1);},x,y_plot);
		plot->write_ave_err(0.0,pr_chir[iout].ave_err());
	      }
	    
	    if(plot) delete plot;
	  }
      }
}

void perens_t::assemble_pr_meslep_QED_greenfunctions()
{
  cout<<"Assembling pr_meslep QED greenfunctions"<<endl;
  
  for(size_t im_in=0;im_in<nm;im_in++)
    for(size_t r_in=0;r_in<nr;r_in++)
      for(size_t im_ou=0;im_ou<nm;im_ou++)
	for(size_t r_ou=0;r_ou<nr;r_ou++)
	  {
	    const size_t im_r_in=im_r_ind({im_in,r_in});
	    const size_t im_r_ou=im_r_ind({im_ou,r_ou});
	    
#pragma omp parallel for
	    for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
	      {
		using namespace meslep;
		
		for(size_t ijack=0;ijack<=njacks;ijack++)
		    for(size_t iop=0;iop<nZop;iop++)
		      for(size_t iproj=0;iproj<nZop;iproj++)
			{
			  using namespace pr_meslep;
			  
			  auto ml=pr_meslep_ins_getter(im_ou,r_ou,im_in,r_in,iop,iproj,imeslepmom);
			  
			  //QED correction
			  switch(pars::use_QED)
			    {
			    case 0:
			      break;
			      /////////////////////////////////////////////////////////////////
			    case 1:
			      ml(QED)[ijack]=
				ml(EX)[ijack]*q_in*q_ou+
				
				(ml(NA_IN)[ijack]*q_in+
				 ml(NA_OU)[ijack]*q_ou)*ql+
				
				(ml(PH_IN)[ijack]+
				 ml(CR_IN)[ijack]*deltam_cr[im_r_in][ijack]+
				 ml(TM_IN)[ijack]*deltam_tm[im_r_in][ijack])*q_in*q_in+
				
				(ml(PH_OU)[ijack] +
				 ml(CR_OU)[ijack]*deltam_cr[im_r_ou][ijack]+
				 ml(TM_OU)[ijack]*deltam_tm[im_r_ou][ijack])*q_ou*q_ou
				;
			      break;
			    case 2:
			      ml(QED)[ijack]=
				ml(EX)[ijack]*q_in*q_ou+
				
				(ml(NA_IN)[ijack]*q_in+
				 ml(NA_OU)[ijack]*q_ou)*ql+
				
				ml(QED_IN)[ijack]*q_in*q_in+
				
				ml(QED_OU)[ijack]*q_ou*q_ou
				;
			      break;
			    }
			}
		  }
	      }
}
