#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_MESLEP
 #include <MOM2/Zmeslep.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

using namespace meslep;

mesloop_t perens_t::build_mesloop(const vector<mom_conf_lprops_t> &props_lep) const
{
  //! projected lepton propagator with insertion
  mesloop_t mesloop(ilistGl_ilistpGl_iclust_ind.max());
  
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
      auto p0=lepGamma[5]*pl.LO.adjoint()*lepGamma[5];
      auto pF=lepGamma[5]*pl.F.adjoint()*lepGamma[5];
      auto pFamp=pF*p0.inverse();
      auto pr=(lepGamma[ipGl]*(lepGamma[0]+psign*lepGamma[5])).adjoint()/2.0;
      mesloop.LO[i]=(op*pr).toDense().trace()/4.0;   //for test: this must be 1 if iGl==ipGl
      mesloop.F[i]=(op*pFamp*pr).trace()/4.0; //normalization for the single gamma
      //cout<<"p0 "<<p0(0,0)<<"\npF "<<pF(0,0)<<"\npFamp "<<pFamp(0,0)<<"\npr "<<((lprop_t)pr)(0,0)<<"\nml "<<mesloop.LO[i]<<"\nmlF "<<mesloop.F[i]<<endl;
    }
  
  return mesloop;
}

void perens_t::build_all_mr_gmeslep_jackkniffed_verts(jmeslep_vert_t &j,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou,
					    const vector<mom_conf_lprops_t> &props_lep) const
{
  const mesloop_t mesloop=build_mesloop(props_lep);
  
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
	 {&j.CR_CT_IN,&mesloop.LO,&p_in.P,&p_ou.LO}, //critical counterterm_in
	 {&j.CR_CT_OU,&mesloop.LO,&p_in.LO,&p_ou.P}, //critical counterterm_ou
	 //
	 //
	 {&j.TM_CT_IN,&mesloop.LO,&p_in.S,&p_ou.LO}, //twisted counterterm_in
	 {&j.TM_CT_OU,&mesloop.LO,&p_in.LO,&p_ou.S}, //twisted counterterm_ou
	 };
       
       const Zop_t &zop=zops[iop];
       
       //create the vertex
       const size_t im_r_im_r_iop_ilistpGl=im_r_im_r_iop_ilistpGl_ind({im_in,r_in,im_ou,r_ou,iop,ilistpGl});
       for(auto &o : list)
	 {
	   const int    sign=zop.Qg5_sign;
	   const int    norm=zop.norm;
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
	       jvert[iclust]+=c*norm;
	       
	       if(1)
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
	   //cout<<" "<<jvert[iclust]<<endl;
	 }
     }
}

djvec_t perens_t::compute_proj_meslep(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou) const
{
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
      
      const double pnorm=12.0*zops[iproj].pnorm*2.0;
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
	  auto contr=(amp_vert*projector).trace().real()/pnorm;
	  out+=contr;
	  
	  if(0)
	  cout<<" amputating, i: "<<i<<", iop: "<<iop<<", ilistpGl: "<<ilistpGl<<", iproj: "<<iproj<<", ijack: "<<ijack<<
	    ", im_r_im_r_iop_ilistpGl(jverts): "<<im_r_im_r_iop_ilistpGl<<", im_r_im_r_iop_iproj(out): "<<im_r_im_r_iop_iproj<<
	    ", prop_inv_in: "<<prop_inv_in(0,0)<<", vert: "<<vert(0,0)<<", prop_inv_ou: "<<prop_inv_ou(0,0)<<", contr: "<<contr<<endl;
	  }
      
      //cout<<"amputated, im_in: "<<im_in<<", r_in: "<<r_in<<", im_ou: "<<im_ou<<", r_ou: "<<r_ou<<", iop: "<<iop<<", iproj: "<<iproj<<", ijack: "<<ijack<<", out: "<<out<<endl;
    }
  
  return pr;
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
      
      vector<jm_r_mom_qprops_t> jprops_in(im_r_ind.max());                //!< jackknived props
      vector<jm_r_mom_qprops_t> jprops_ou(im_r_ind.max());                //!< jackknived props
      jmeslep_vert_t jmeslep_verts(im_r_im_r_ilistGl_ipGl_ind.max());     //!< jackknived vertex
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
	    vector<m_r_mom_conf_qprops_t> props_in=read_all_qprops_mom(qfiles,i_in_clust_ihit,mom_in);
	    vector<m_r_mom_conf_qprops_t> props_ou=(read_ou?read_all_qprops_mom(qfiles,i_in_clust_ihit,mom_ou):props_in);
	    vector<mom_conf_lprops_t> props_lep=read_all_lprops_mom(lfiles,i_in_clust_ihit,momlep);
	    read_time.stop();
	    
	    //put the charges in place
	    using namespace meslep;
	    incorporate_charge(props_in,q_in);
	    incorporate_charge(props_ou,q_ou);
	    incorporate_charge(props_lep,ql);
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops_in,props_in);
	    build_all_mr_jackkniffed_qprops(jprops_ou,props_ou);
	    build_props_time.stop();
	    
	    //build all meslep_verts
	    build_meslep_verts_time.start();
	    build_all_mr_gmeslep_jackkniffed_verts(jmeslep_verts,props_in,props_ou,props_lep);
	    build_meslep_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops_in);
      clusterize_all_mr_jackkniffed_qprops(jprops_ou);
      jmeslep_verts.clusterize_all(clust_size,im_r_im_r_iop_ilistpGl_ind,deltam_cr,deltam_tm);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv_in;      //!< inverse propagator_in
      vector<jqprop_t> jprop_inv_ou;      //!< inverse propagator_ou
      vector<jqprop_t> jprop_QED_inv_in;  //!< inverse propagator_in with em insertion
      vector<jqprop_t> jprop_QED_inv_ou;  //!< inverse propagator_ou with em insertion
      
      get_inverse_propagators(jprop_inv_in,jprop_QED_inv_in,jprops_in);
      get_inverse_propagators(jprop_inv_ou,jprop_QED_inv_ou,jprops_ou);
      
      proj_time.start();
      
      const auto &j=jmeslep_verts;
      
      djvec_t pr_LO,pr_QED_amp_QCD,pr_QCD_amp_QED_in,pr_QCD_amp_QED_ou;
      for(auto &p : vector<tuple<djvec_t*,const vector<jqprop_t>*,const vector<jqprop_t>*,const vector<jqprop_t>*>>{
	  {&pr_LO,              &jprop_inv_in,     &j.LO,  &jprop_inv_ou},
	  {&pr_QED_amp_QCD,     &jprop_inv_in,     &j.QED, &jprop_inv_ou},
	  {&pr_QCD_amp_QED_in,  &jprop_QED_inv_in, &j.LO,  &jprop_inv_ou},
	  {&pr_QCD_amp_QED_ou,  &jprop_inv_in,     &j.LO,  &jprop_QED_inv_ou}
	    })
	{
	  auto &out=*get<0>(p);
	  auto &p_in=*get<1>(p);
	  auto &v=*get<2>(p);
	  auto &p_ou=*get<3>(p);
	  out=compute_proj_meslep(p_in,v,p_ou);
	}
      proj_time.stop();
      
      //! an index running on all packed combo, and momenta
      const index_t all_imeslepmom_ind({{"All",im_r_im_r_iop_iproj_ind.max()},{"meslepmom",meslepmoms().size()}});
      
      //store
      for(size_t iall=0;iall<im_r_im_r_iop_iproj_ind.max();iall++)
	{
	  pr_meslep[all_imeslepmom_ind({iall,imeslepmom})]=
	    +pr_LO[iall];
	  pr_meslep_QED[all_imeslepmom_ind({iall,imeslepmom})]=
	    +pr_QED_amp_QCD[iall]
	    -pr_QCD_amp_QED_in[iall]
	    -pr_QCD_amp_QED_ou[iall]
	    ;
	}
    }
}

vector<perens_t::task_t> perens_t::get_pr_meslep_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_pr_meslep,in_pr_meslep_QED;
  for(auto &e : ens)
    {
      in_pr_meslep.push_back(&e->pr_meslep);
      if(pars::use_QED) in_pr_meslep_QED.push_back(&e->pr_meslep_QED);
    }
  vector<task_t> pr_meslep_tasks={{&pr_meslep,in_pr_meslep,im_r_im_r_iop_iproj_imeslepmom_ind,"pr_meslep",QCD_task}};
  if(pars::use_QED) pr_meslep_tasks.push_back({&pr_meslep_QED,in_pr_meslep_QED,im_r_im_r_iop_iproj_imeslepmom_ind,"pr_meslep_QED",QED_task});
  
  return pr_meslep_tasks;
}

vector<perens_t::task_t> perens_t::get_Zmeslep_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_Zmeslep,in_Zmeslep_QED_rel;
  for(auto &e : ens)
    {
      in_Zmeslep.push_back(&e->Zmeslep);
      if(pars::use_QED) in_Zmeslep_QED_rel.push_back(&e->Zmeslep_QED_rel);
    }
  vector<task_t> Zmeslep_tasks={{&Zmeslep,in_Zmeslep,im_r_im_r_iop_iproj_imeslepmom_ind,"Zmeslep",QCD_task}};
  if(pars::use_QED) Zmeslep_tasks.push_back({&Zmeslep_QED_rel,in_Zmeslep_QED_rel,im_r_im_r_iop_iproj_imeslepmom_ind,"Zmeslep_QED_rel",QED_task});
  
  return Zmeslep_tasks;
}

void perens_t::compute_Zmeslep()
{
  for(size_t im_in=0;im_in<nm;im_in++)
    for(size_t r_in=0;r_in<nr;r_in++)
      for(size_t im_ou=0;im_ou<nm;im_ou++)
	for(size_t r_ou=0;r_ou<nr;r_ou++)
#pragma omp parallel for
	  for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
	    {
	      using namespace meslep;
	      
	      const size_t ilinmom_in=meslepmoms()[imeslepmom][1];
	      const size_t ilinmom_ou=meslepmoms()[imeslepmom][2];
	      const size_t im_r_in_ilinmom_in=im_r_ilinmom_ind({im_in,r_in,ilinmom_in});
	      const size_t im_r_ou_ilinmom_ou=im_r_ilinmom_ind({im_ou,r_ou,ilinmom_ou});
	      
	      for(size_t ijack=0;ijack<=njacks;ijack++)
		{
		  using Zmeslep_t=Matrix<double,nZop,nZop>;
		  Zmeslep_t Gamma_meslep_combo;
		  Zmeslep_t Gamma_QED_meslep_combo;
		  
		  for(size_t iop=0;iop<nZop;iop++)
		    for(size_t iproj=0;iproj<nZop;iproj++)
		      {
			const size_t im_r_im_r_iop_iproj_imeslepmom=im_r_im_r_iop_iproj_imeslepmom_ind({im_in,r_in,im_ou,r_ou,iop,iproj,imeslepmom});
			
			//Gamma LO and correction
			Gamma_meslep_combo(iop,iproj)=pr_meslep[im_r_im_r_iop_iproj_imeslepmom][ijack];
			if(pars::use_QED) Gamma_QED_meslep_combo(iop,iproj)=pr_meslep_QED[im_r_im_r_iop_iproj_imeslepmom][ijack];
		      }
		  
		  const Zmeslep_t Gamma_meslep_combo_inv=Gamma_meslep_combo.inverse();
		  
		  auto Zq_contr=sqrt(Zq_sig1[im_r_in_ilinmom_in][ijack]*Zq_sig1[im_r_ou_ilinmom_ou][ijack]);
		  auto Zq_QED_rel_contr=
		    0.5*
		    (Zq_sig1_QED[im_r_in_ilinmom_in][ijack]/Zq_sig1[im_r_in_ilinmom_in][ijack]*sqr(meslep::q_in)+
		     Zq_sig1_QED[im_r_ou_ilinmom_ou][ijack]/Zq_sig1[im_r_ou_ilinmom_ou][ijack]*sqr(meslep::q_ou));
		  
		  auto Z_LO=Zq_contr*Gamma_meslep_combo_inv;
		  auto Z_QED_rel=(Zq_QED_rel_contr*Zmeslep_t::Identity()
				  -Gamma_QED_meslep_combo*Gamma_meslep_combo_inv);
		  
		  for(size_t iop=0;iop<nZop;iop++)
		    for(size_t iproj=0;iproj<nZop;iproj++)
		      {
			const size_t im_r_im_r_iop_iproj_imeslepmom=im_r_im_r_iop_iproj_imeslepmom_ind({im_in,r_in,im_ou,r_ou,iop,iproj,imeslepmom});
			
			Zmeslep[im_r_im_r_iop_iproj_imeslepmom][ijack]=Z_LO(iop,iproj);
			Zmeslep_QED_rel[im_r_im_r_iop_iproj_imeslepmom][ijack]=Z_QED_rel(iop,iproj);
		      }
		}
	    }
}

void perens_t::plot_Zmeslep(const string &suffix)
{
  for(const auto &t : this->get_Zmeslep_tasks())
      {
	//decript tuple
	const djvec_t &Z=*t.out;
	const string &tag=t.tag;
	
	for(size_t iop=0;iop<nZbil;iop++)
	  for(size_t iproj=0;iproj<nZbil;iproj++)
	    {
	      grace_file_t out(dir_path+"/plots/"+tag+"_"+to_string(iop)+"_"+to_string(iproj)+(suffix!=""?("_"+suffix):string(""))+".xmg");
	      
	      //write mass by mass, only half of the combos
	      for(size_t im_in=0;im_in<nm;im_in++)
		for(size_t im_ou=im_in;im_ou<nm;im_ou++)
		  for(size_t r=0;r<nr;r++)
		    {
		      out.new_data_set();
		      
		      for(size_t imom=0;imom<meslepmoms().size();imom++)
			{
			  const double p2tilde=all_moms[meslepmoms()[imom][0]].p(L).tilde().norm2();
			  size_t i=im_r_im_r_iop_iproj_imeslepmom_ind({im_in,r,im_ou,r,iop,iproj,imom});
			  out.write_ave_err(p2tilde,Z[i].ave_err());
			}
		    }
	    }
      }
}

void perens_t::average_r_Zmeslep(perens_t &out) const
{
  cout<<"Averaging r for meslep"<<endl;
  
  for(auto &t : out.get_meslep_tasks({this}))
    {
      cout<<" "<<t.tag<<endl;
      
      const djvec_t &pr=*t.in.front();
      djvec_t &pr_rave=*t.out;
      
      for(size_t out_i=0;out_i<out.im_r_im_r_iop_iproj_imeslepmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_iop_iproj_imeslepmom_comp=out.im_r_im_r_iop_iproj_imeslepmom_ind(out_i);
	  vector<size_t> im_r_im_r_iop_iproj_imeslepmom_comp=out_im_r_im_r_iop_iproj_imeslepmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_im_r_iop_iproj_imeslepmom_comp[1]=
		im_r_im_r_iop_iproj_imeslepmom_comp[3]=r;
	      const size_t i=im_r_im_r_iop_iproj_imeslepmom_ind(im_r_im_r_iop_iproj_imeslepmom_comp);
	      pr_rave[out_i]+=pr[i];
	    }
	  pr_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_Zmeslep(perens_t &out,const vector<vector<size_t>> &equiv_meslepmom_combos) const
{
  for(size_t i=0;i<out.im_r_im_r_iop_iproj_imeslepmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_iop_iproj_imeslepmom_comp=out.im_r_im_r_iop_iproj_imeslepmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_iop_iproj_imeslepmom_comp[6];
      
      for(const auto &t : out.get_meslep_tasks({this}))
	{
	  djack_t &ave=(*t.out)[i];
	  ave=0.0;
	  for(const size_t ieq_mom : equiv_meslepmom_combos[out_imom_combo])
	    {
	      vector<size_t> im_r_im_r_iop_iproj_imeslepmom_comp=out_im_r_im_r_iop_iproj_imeslepmom_comp;
	      im_r_im_r_iop_iproj_imeslepmom_comp[6]=ieq_mom;
	      const size_t ieq=im_r_im_r_iop_iproj_imeslepmom_ind(im_r_im_r_iop_iproj_imeslepmom_comp);
	      ave+=(*t.in.front())[ieq];
	    }
	  ave/=equiv_meslepmom_combos[out_imom_combo].size();
	}
    }
}

void perens_t::val_chir_extrap_Zmeslep(perens_t &out) const
{
  for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
    for(auto &t : out.get_meslep_tasks({this}))
      for(size_t iop=0;iop<meslep::nZop;iop++)
	for(size_t iproj=0;iproj<meslep::nZop;iproj++)
	  {
	    const djvec_t &pr=*t.in.front();
	    djvec_t &pr_chir=*t.out;
	    const string &tag=t.tag;
	    
	    //check if we need to subtract the pole
	    const bool sub_pole=(iop==2 or iop==3 or iproj==2 or iproj==3);
	    const size_t coeff_to_take=(sub_pole?1:0);
	    
	    //open the plot file if needed
	    const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+to_string(iop)+"_"+to_string(iproj)+"_meslepmom_"+to_string(imeslepmom)+".xmg";
	    grace_file_t *plot=nullptr;
	    if(imeslepmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	    
	    for(size_t r1=0;r1<nr;r1++)
	      for(size_t r2=0;r2<nr;r2++)
		{
		  //slice m and fit it
		  djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
		  vector<double> x(nm*(nm+1)/2);
		  int i=0;
		  for(size_t im1=0;im1<nm;im1++)
		    for(size_t im2=im1;im2<nm;im2++)
		      {
			//compute mass sum
			if(pars::chir_extr_method==chir_extr::MQUARK) x[i]=am[im1]+am[im2];
			else                                          x[i]=sqr(meson_mass[im_im_ind({im1,im2})].ave());
			//compute y and y_plot
			y_plot[i]=pr[im_r_im_r_iop_iproj_imeslepmom_ind({im1,r1,im2,r2,iop,iproj,imeslepmom})];
			//fit x*y if pole present
			if(sub_pole) y[i]=x[i]*y_plot[i];
			else         y[i]=y_plot[i];
			//increment the number of mass combos
			i++;
		      }
		  
		  //fit, store and write the result
		  djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1));
		  
		  if(std::isnan(coeffs[0][0])) coeffs=0.0;
		  
		  const size_t iout=out.im_r_im_r_iop_iproj_imeslepmom_ind({0,r1,0,r2,iop,iproj,imeslepmom});
		  pr_chir[iout]=coeffs[coeff_to_take];
		  if(plot!=nullptr)
		    {
		      auto xminmax=minmax_element(x.begin(),x.end());
		      double xmin=*xminmax.first*0.9;
		      double xmax=*xminmax.second*1.1;
		      write_fit_plot(*plot,xmin,xmax,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
		      plot->write_ave_err(0.0,pr_chir[iout].ave_err());
		    }
		}
	    if(plot) delete plot;
	  }
}
