#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_INGREDIENTS
 #include <ingredients.hpp>

#include <deltam_cr.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

void ingredients_t::set_pars_for_scratch()
{
  using namespace reno_scheme;
  
  _nm=glb::nm;
  _nr=glb::nr;
  _am=glb::am;
  
  switch(scheme)
    {
    case RI_MOM:
      set_ri_mom_moms();
      break;
    case SMOM:
      set_smom_moms();
      break;
    }
}

void ingredients_t::set_ri_mom_moms()
{
  for(size_t imom=0;imom<glb_moms.size();imom++)
    if(filt_moms[imom])
      {
	const size_t linmom=linmoms.size();
	linmoms.push_back({imom});
	bilmoms.push_back({imom,linmom,linmom});
    }
}

void ingredients_t::set_smom_moms()
{
  CRASH("Fix: list linmom only useful for smom, set other two components of bimomd relatively to linmom, not glb");
  
  for(size_t imom=0;imom<glb_moms.size();imom++)
    if(filt_moms[imom])
      linmoms.push_back({imom});
  
  const double tol=1e-10;
  for(size_t i=0;i<glb_moms.size();i++)
    if(filt_moms[i])
      {
	//get norm of p[i]
	p_t pi=glb_moms[i].p(L);
	double pi2=pi.norm2();
	
	for(size_t j=0;j<glb_moms.size();j++)
	  if(filt_moms[j])
	    {
	      //get norm of p[j]
	      p_t pj=glb_moms[j].p(L);
	      double pj2=pj.norm2();
	      
	      //check that norm of the two incoming vector is the same
	      if(2.0*fabs(pi2-pj2)<(pi2+pj2)*tol)
		{
		  //sum and get norm
		  imom_t momk;
		  for(size_t mu=0;mu<NDIM;mu++)
		    momk[mu]=glb_moms[i][mu]+glb_moms[j][mu];
		  double pk2=momk.p(L).norm2();
		  
		  //debug info
		  //cerr<<"pi2: "<<pi2<<" pk2: "<<pk2<<" "<<2.0*fabs(pi2-pk2)/(pi2+pk2)<<" "<<pi2<<" "<<pj2<<" "<<pk2<<endl;
		  
		  if(2.0*fabs(pi2-pk2)<(pi2+pk2)*tol)
		    {
		      //search in list
		      auto posk=find(glb_moms.begin(),glb_moms.end(),momk);
		      //take position if found
		      if(posk!=glb_moms.end())
			{
			  const size_t k=distance(glb_moms.begin(),posk);
			  //inform and add to the list
			  //cout<<"Found smom pair: "<<i<<" "<<glb_moms[i]<<" + j "<<glb_moms[j]<<" = "<<k<<" "<<momk<<endl;
			  bilmoms.push_back({k,i,j});
			}
		      // else
		      //    cout<<"Unable to find it"<<momk<<"="<<glb_moms[i]<<"+"<<glb_moms[j]<<endl;
		    }
		}
	    }
      }
  
  cout<<"Number of smom pairs: "<<bilmoms.size()<<endl;
}

void ingredients_t::mom_compute_prop()
{
  vector<raw_file_t> files=setup_read_all_props_mom(conf_list);
  
  for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
    {
      vector<jm_r_mom_props_t> jprops(glb::im_r_ind.max()); //!< jackknived props
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", hit "<<ihit+1<<"/"<<nhits<<", momentum "<<ilinmom+1<<"/"<<linmoms.size()<<endl;
	    read_time.start();
	    const vector<m_r_mom_conf_props_t> props=read_all_props_mom(files,i_in_clust_ihit,linmoms[ilinmom][0]);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackknifed_props(jprops,props,use_QED,glb::im_r_ind,deltam_cr);
	    build_props_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      clust_time.stop();
      
      vector<jprop_t> jprop_inv(glb::im_r_ind.max()); //!< inverse propagator
      vector<jprop_t> jprop_EM_inv(glb::im_r_ind.max()); //!< inverse propagator with em insertion
#pragma omp parallel for reduction(+:invert_time,Zq_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=glb::im_r_ind({im,r});
	  const size_t im_r_ilinmom=im_r_ilinmom_ind({im,r,ilinmom});
	  
	  //compute inverse
	  invert_time.start();
	  const prop_t prop_inv=jprops[im_r].LO[ijack].inverse();
	  //if(im_r_ijack==0) cout<<jprops[im_r].LO[ijack](0,0)<<endl;
	  jprop_inv[im_r][ijack]=prop_inv;
	  invert_time.stop();
	  
	  //compute Zq
	  Zq_time.start();
	  Zq[im_r_ilinmom][ijack]=compute_Zq(prop_inv,ilinmom);
	  Zq_sig1[im_r_ilinmom][ijack]=compute_Zq_sig1(prop_inv,ilinmom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(use_QED)
	    {
	      invert_time.start();
	      prop_t prop_EM_inv=prop_inv*jprops[im_r].EM[ijack]*prop_inv;
	      jprop_EM_inv[im_r][ijack]=prop_EM_inv;
	      invert_time.stop();
	      
	      Zq_time.start();
	      Zq_sig1_EM[im_r_ilinmom][ijack]=-compute_Zq_sig1(prop_EM_inv,ilinmom);
	      Zq_time.stop();
	    }
	}
    }
}

void ingredients_t::mom_compute_bil()
{
  vector<raw_file_t> files=setup_read_all_props_mom(conf_list);
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    {
      const size_t imom1=bilmoms[ibilmom][1];
      const size_t imom2=bilmoms[ibilmom][2];
      const bool read2=(imom1!=imom2);
      
      vector<jm_r_mom_props_t> jprops1(glb::im_r_ind.max()); //!< jackknived props
      vector<jm_r_mom_props_t> jprops2(glb::im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED);  //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", hit "<<ihit+1<<"/"<<nhits<<", momentum combo "<<ibilmom+1<<"/"<<bilmoms.size()<<endl;
	    read_time.start();
	    const vector<m_r_mom_conf_props_t> props1=read_all_props_mom(files,i_in_clust_ihit,imom1);
	    const vector<m_r_mom_conf_props_t> props2=(read2)?read_all_props_mom(files,i_in_clust_ihit,imom2):props1;
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackknifed_props(jprops1,props1,use_QED,glb::im_r_ind,deltam_cr);
	    build_all_mr_jackknifed_props(jprops2,props2,use_QED,glb::im_r_ind,deltam_cr);
	    build_props_time.stop();
	    
	    build_verts_time.start();
	    build_all_mr_gbil_jackknifed_verts(jverts,props1,props2,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED,deltam_cr);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops1,use_QED,clust_size);
      clusterize_all_mr_jackknifed_props(jprops2,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
      clust_time.stop();
      
      vector<jprop_t> jprop_inv1(glb::im_r_ind.max()); //!< inverse propagator1
      vector<jprop_t> jprop_inv2(glb::im_r_ind.max()); //!< inverse propagator2
      vector<jprop_t> jprop_EM_inv1(glb::im_r_ind.max()); //!< inverse propagator1 with em insertion
      vector<jprop_t> jprop_EM_inv2(glb::im_r_ind.max()); //!< inverse propagator2 with em insertion
#pragma omp parallel for reduction(+:invert_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=glb::im_r_ind({im,r});
	  
	  //compute inverse
	  invert_time.start();
	  prop_t prop_inv1=jprop_inv1[im_r][ijack]=jprops1[im_r].LO[ijack].inverse();
	  prop_t prop_inv2=jprop_inv2[im_r][ijack]=jprops2[im_r].LO[ijack].inverse();
	  //if(im_r_ijack==0) cout<<jprops[im_r].LO[ijack](0,0)<<endl;
	  invert_time.stop();
	  
	  //do the same with QED
	  if(use_QED)
	    {
	      invert_time.start();
	      jprop_EM_inv1[im_r][ijack]=prop_inv1*jprops1[im_r].EM[ijack]*prop_inv1;
	      jprop_EM_inv2[im_r][ijack]=prop_inv2*jprops2[im_r].EM[ijack]*prop_inv2;
	      invert_time.stop();
	    }
	}
      
      //! an index running on all packed combo, and momenta
      const index_t all_ibilmom_ind({{"All",im_r_im_r_iZbil_ind.max()},{"bilmom",bilmoms.size()}});
      
      proj_time.start();
      djvec_t pr_bil_temp=compute_proj_bil(jprop_inv1,jverts.LO,jprop_inv2,glb::im_r_ind);
      
      //QED
      djvec_t pr_bil_QED_temp;
      if(use_QED)
	{
	  const djvec_t pr_bil_EM=compute_proj_bil(jprop_inv1,jverts.EM,jprop_inv2,glb::im_r_ind);
	  const djvec_t pr_bil_a=compute_proj_bil(jprop_EM_inv1,jverts.LO,jprop_inv2,glb::im_r_ind);
	  const djvec_t pr_bil_b=compute_proj_bil(jprop_inv1,jverts.LO,jprop_EM_inv2,glb::im_r_ind);
	  pr_bil_QED_temp=pr_bil_a+pr_bil_b-pr_bil_EM;
	}
      
      //store
      vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil,&pr_bil_temp}};
      if(use_QED) store_task.push_back({&pr_bil_QED,&pr_bil_QED_temp});
      for(auto &t : store_task)
	for(size_t iall=0;iall<im_r_im_r_iZbil_ind.max();iall++)
	  (*t.first)[all_ibilmom_ind({iall,ibilmom})]=(*t.second)[iall];
      
      proj_time.stop();
    }
}

void ingredients_t::ri_mom()
{
  mom_compute_prop();
  mom_compute_bil();
}

void ingredients_t::smom()
{
  ri_mom();
}

void ingredients_t::bin_read(raw_file_t &file)
{
  for(djvec_t *o: {&Zq,&Zq_sig1,&Zq_sig1_EM,&pr_bil,&pr_bil_QED}) o->bin_read(file);
}

void ingredients_t::bin_write(raw_file_t &file) const
{
  for(const djvec_t *o: {&Zq,&Zq_sig1,&Zq_sig1_EM,&pr_bil,&pr_bil_QED}) o->bin_write(file);
}

void ingredients_t::set_indices()
{
  im_r_ind.set_ranges({{"m",_nm},{"r",_nr}});
  im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"iZbil",nZbil}});;
  r_ilinmom_ind.set_ranges({{"r",_nr},{"linmom",linmoms.size()}});
  im_r_ilinmom_ind.set_ranges({{"m",_nm},{"r",_nr},{"linmom",linmoms.size()}});
  iZbil_ibilmom_ind.set_ranges({{"iZbil",nZbil},{"bilmom",bilmoms.size()}});
  im_r_im_r_iZbil_ibilmom_ind=im_r_im_r_iZbil_ind*index_t({{"bilmom",bilmoms.size()}});
}

void ingredients_t::allocate()
{
  Zq.resize(im_r_ilinmom_ind.max());
  Zq_sig1.resize(im_r_ilinmom_ind.max());
  if(use_QED) Zq_sig1_EM.resize(im_r_ilinmom_ind.max());
  
  pr_bil.resize(im_r_im_r_iZbil_ibilmom_ind.max());
  if(use_QED) pr_bil_QED.resize(im_r_im_r_iZbil_ibilmom_ind.max());
  
  Zbil.resize(im_r_im_r_iZbil_ibilmom_ind.max());
  if(use_QED) Zbil_QED.resize(im_r_im_r_iZbil_ibilmom_ind.max());
}

void ingredients_t::create_from_scratch(const string ingredients_path)
{
  set_pars_for_scratch();
  
  set_indices();
  
  allocate();
  
  //if ingredients exists read it, otherwise compute it
  if(file_exists(ingredients_path)) bin_read(ingredients_path);
  else
    {
      prepare_list_of_confs();
      
      using namespace reno_scheme;
      
      switch(scheme)
	{
	case RI_MOM:
	  ri_mom();
	  break;
	case SMOM:
	  smom();
	  break;
	}
      bin_write(ingredients_path);
    }
  
  compute_Zbil();
}

ingredients_t ingredients_t::average_r(const bool recompute_Zbil) const
{
  ingredients_t out;
  
  out._nr=1;
  out._nm=_nm;
  out._am=_am;
  out.linmoms=linmoms;
  out.bilmoms=bilmoms;
  out.Zbil_computed=Zbil_computed;
  
  out.set_indices();
  out.allocate();
  
  //average Zq
  for(auto &t : get_Zq_tasks(out))
    {
      const djvec_t &Zq=*t.in;
      djvec_t &Zq_rave=*t.out;
      
      for(size_t out_i=0;out_i<out.im_r_ilinmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(out_i);
	  vector<size_t> im_r_ilinmom_comp=out_im_r_ilinmom_comp;
	  
	  Zq_rave[out_i]=0.0;
	  for(size_t r=0;r<_nr;r++)
	    {
	      im_r_ilinmom_comp[1]=r;
	      const size_t i=im_r_ilinmom_ind(im_r_ilinmom_comp);
	      Zq_rave[out_i]+=Zq[i];
	    }
	  Zq_rave[out_i]/=_nr;
	}
    }
  
  //average pr_bil
  for(auto &t : get_pr_bil_tasks(out))
    {
      const djvec_t &pr=*t.in;
      djvec_t &pr_rave=*t.out;
      
      for(size_t out_i=0;out_i<out.im_r_im_r_iZbil_ibilmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_iZbil_ibilmom_comp=out.im_r_im_r_iZbil_ibilmom_ind(out_i);
	  vector<size_t> im_r_im_r_iZbil_ibilmom_comp=out_im_r_im_r_iZbil_ibilmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<_nr;r++)
	    {
	      im_r_im_r_iZbil_ibilmom_comp[1]=
		im_r_im_r_iZbil_ibilmom_comp[3]=r;
	      const size_t i=im_r_im_r_iZbil_ibilmom_ind(im_r_im_r_iZbil_ibilmom_comp);
	      pr_rave[out_i]+=pr[i];
	    }
	  pr_rave[out_i]/=_nr;
	}
    }
   
  return out;
}

void ingredients_t::compute_Zbil()
{
  Zbil_computed=true;
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(size_t im_r_im_r_iZbil=0;im_r_im_r_iZbil<im_r_im_r_iZbil_ind.max();im_r_im_r_iZbil++)
      {
	const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
	const vector<size_t> im_r1_comp=subset(im_r_im_r_iZbil_comp,0,2);
	const vector<size_t> im_r2_comp=subset(im_r_im_r_iZbil_comp,2,4);
	const size_t iZbil=im_r_im_r_iZbil_comp[4];
	const size_t ilinmom1=bilmoms[ibilmom][1];
	const size_t ilinmom2=bilmoms[ibilmom][2];
	const size_t im_r1_ilinmom1=im_r_ilinmom_ind(concat(im_r1_comp,ilinmom1));
	const size_t im_r2_ilinmom2=im_r_ilinmom_ind(concat(im_r2_comp,ilinmom2));
	
	const size_t im_r_im_r_iZbil_ibilmom=im_r_im_r_iZbil_ibilmom_ind(concat(im_r1_comp,im_r2_comp,vector<size_t>({iZbil,ibilmom})));
	
	Zbil[im_r_im_r_iZbil_ibilmom]=
	  sqrt(Zq_sig1[im_r1_ilinmom1]*Zq_sig1[im_r2_ilinmom2])/pr_bil[im_r_im_r_iZbil_ibilmom];
	
	if(use_QED)
	  {
	    Zbil_QED[im_r_im_r_iZbil_ibilmom]=
	      pr_bil_QED[im_r_im_r_iZbil_ibilmom]/pr_bil[im_r_im_r_iZbil_ibilmom]+
	      (Zq_sig1_EM[im_r1_ilinmom1]/Zq_sig1[im_r1_ilinmom1]+Zq_sig1_EM[im_r2_ilinmom2]/
	       Zq_sig1[im_r2_ilinmom2])/2.0;
	  }
      }
}

ingredients_t ingredients_t::chir_extrap() const
{
  ingredients_t out;
  out._nr=_nr;
  out._nm=1;
  out._am={0.0};
  out.linmoms=linmoms;
  out.bilmoms=bilmoms;
  out.Zbil_computed=Zbil_computed;
  
  out.set_indices();
  out.allocate();
  
  //extrapolate to chiral limit Zq
  for(auto &t : get_Zq_tasks(out))
    {
      const djvec_t &Zq=*t.in;
      djvec_t &Zq_chir=*t.out;
      const string &tag=t.tag;
      
      for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
	{
	  //open the plot file if needed
	  const string plot_path="plots/chir_extr_"+tag+"_mom_"+to_string(ilinmom)+".xmg";
	  grace_file_t *plot=nullptr;
	  if(ilinmom%print_each_mom==0) plot=new grace_file_t(plot_path);
	  
	  for(size_t r=0;r<_nr;r++)
	    {
	      //slice m
	      djvec_t y(_nm);
	      for(size_t im=0;im<_nm;im++) y[im]=Zq[im_r_ilinmom_ind({im,r,ilinmom})];
	      
	      //fit, store and write the result
	      djvec_t coeffs=poly_fit(_am,y,1,am_min,am_max);
	      Zq_chir[out.im_r_ilinmom_ind({0,r,ilinmom})]=coeffs[0];
	      if(plot!=nullptr)
		{
		  write_fit_plot(*plot,0,am_max,bind(poly_eval<djvec_t>,coeffs,_1),_am,y);
		  plot->write_ave_err(0,coeffs[0].ave_err());
		}
	    }
	  
	  if(plot!=nullptr) delete plot;
	}
    }
  
  //extrapolate to chiral limit bil
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(auto &t : get_pr_bil_tasks(out))
	for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	  {
	    const djvec_t &pr=*t.in;
	    djvec_t &pr_chir=*t.out;
	    const string &tag=t.tag;
	    
	    //check if we need to subtract the pole
	    const bool sub_pole=(iZbil==iZS or iZbil==iZP);
	    const size_t coeff_to_take=(sub_pole?1:0);
	    
	    //open the plot file if needed
	    const string plot_path="plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_bilmom_"+to_string(ibilmom)+".xmg";
	    grace_file_t *plot=nullptr;
	    if(ibilmom%print_each_mom==0) plot=new grace_file_t(plot_path);
	    
	    for(size_t r1=0;r1<_nr;r1++)
	      for(size_t r2=0;r2<_nr;r2++)
		{
		  //slice m and fit it
		  djvec_t y(_nm*(_nm+1)/2),y_plot(_nm*(_nm+1)/2);
		  vector<double> x(_nm*(_nm+1)/2);
		  int i=0;
		  for(size_t im1=0;im1<_nm;im1++)
		    for(size_t im2=im1;im2<_nm;im2++)
		      {
			//compute mass sum
			x[i]=_am[im1]+_am[im2];
			//compute y and y_plot
			y_plot[i]=pr[im_r_im_r_iZbil_ibilmom_ind({im1,r1,im2,r2,iZbil,ibilmom})];
			//fit x*y if pole present
			if(sub_pole) y[i]=x[i]*y_plot[i];
			else         y[i]=y_plot[i];
			//increment the number of mass combos
			i++;
		      }
		  
		  //fit, store and write the result
		  const djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1),2.0*am_min,2.0*am_max);
		  const size_t iout=out.im_r_im_r_iZbil_ibilmom_ind({0,r1,0,r2,iZbil,ibilmom});
		  pr_chir[iout]=coeffs[coeff_to_take];
		  if(plot!=nullptr)
		    {
		      write_fit_plot(*plot,2*am_min,2*am_max,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
		      plot->write_ave_err(0.0,pr_chir[iout].ave_err());
		    }
		}
	    if(plot) delete plot;
	  }
  
  /////////////////////////////////////////////
  
  out.compute_Zbil();
  
  return out;
}

ingredients_t ingredients_t::subtract_Oa2(const bool recompute_Zbil) const
{
  ingredients_t out=(*this);
  
  //subtract Zq
  for(size_t i=0;i<im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> im_r_ilinmom_comp=im_r_ilinmom_ind(i);
      const size_t ilinmom=im_r_ilinmom_comp[2];
      const size_t imom=linmoms[ilinmom][0];
      const imom_t mom=glb_moms[imom];
      const double sub=g2tilde*sig1_a2(act,gf::LANDAU,group::SU3,mom,L);
      out.Zq[i]=Zq[i]-sub;
      out.Zq_sig1[i]=Zq_sig1[ilinmom]-sub;
      if(use_QED)
	{
	  double sub_EM=-1.0*/*coupling gets -1 due to definition of expansion*/sig1_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L);
	  out.Zq_sig1_EM[i]=Zq_sig1_EM[i]-sub_EM*Zq_sig1[i]; //factorization hypotesis
	}
    }
  
  //subtract bil
  for(size_t i=0;i<im_r_im_r_iZbil_ibilmom_ind.max();i++)
    {
      const vector<size_t> im_r_im_r_iZbil_ibilmom_comp=im_r_im_r_iZbil_ibilmom_ind(i);
      const size_t iZbil=im_r_im_r_iZbil_ibilmom_comp[4];
      const size_t ibilmom=im_r_im_r_iZbil_ibilmom_comp[5];
      const size_t imom=bilmoms[ibilmom][0];
      
      imom_t mom=glb_moms[imom];
      double sub;
      sub=g2tilde*pr_bil_a2(act,gf::LANDAU,group::SU3,mom,L,iZbil);
      out.pr_bil[i]=pr_bil[i]-sub;
      if(use_QED)
	{
	  djack_t sub_QED=1.0*pr_bil_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L,iZbil)*pr_bil[i]; //factorization hypotesis
	  out.pr_bil[i]=pr_bil[i]-sub_QED;
	}
    }
   
  if(not Zbil_computed or recompute_Zbil) out.compute_Zbil();
  
  return out;
}

ingredients_t ingredients_t::evolve() const
{
  ingredients_t out=(*this);
  
  //evolve Zq
  for(size_t i=0;i<im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> im_r_ilinmom_comp=im_r_ilinmom_ind(i);
      const size_t ilinmom=im_r_ilinmom_comp[2];
      const size_t imom=linmoms[ilinmom][0];
      
      imom_t mom=glb_moms[imom];
      double p2=mom.p(L).norm2();
      double evolver_Zq=evolution_Zq_to_RIp(Nf,evo_ord,ainv,p2);
      //cout<<"EvolverZq["<<p2<<"]="<<evolver_Zq<<endl;
      
      out.Zq[i]=Zq[i]/evolver_Zq;
      out.Zq_sig1[i]=Zq_sig1[ilinmom]/evolver_Zq;
      if(use_QED) out.Zq_sig1_EM[i]=Zq_sig1_EM[i]/evolver_Zq; //neglecting QED evolution
    }
  
  if(not Zbil_computed) CRASH("Cannot evolve before computing Zbil");
  
  //evolve bil
  for(size_t i=0;i<im_r_im_r_iZbil_ibilmom_ind.max();i++)
    {
      const vector<size_t> im_r_im_r_iZbil_ibilmom_comp=im_r_im_r_iZbil_ibilmom_ind(i);
      const size_t iZbil=im_r_im_r_iZbil_ibilmom_comp[4];
      const size_t ibilmom=im_r_im_r_iZbil_ibilmom_comp[5];
      const size_t imom=bilmoms[ibilmom][0];
      
      const imom_t mom=glb_moms[imom];
      const double p2=mom.p(L).norm2();
      const double evolver_Zbil=evolution_Zbil_to_RIp(iZbil_t_list[iZbil],Nf,evo_ord,ainv,p2);
      //cout<<"EvolverZ"<<Zbil_tag[iZbil]<<"["<<p2<<"]="<<evolver_Zbil<<endl;
      
      out.Zbil[i]=Zbil[i]/evolver_Zbil;
      if(use_QED) out.Zbil[i]=Zbil[i]/evolver_Zbil;
    }
   
  return out;
}

ingredients_t ingredients_t::average_equiv_momenta(const bool recompute_Zbil) const
{
  ingredients_t out;
  out._nm=_nm;
  out._nr=_nr;
  out._am=_am;
  out.Zbil_computed=Zbil_computed;
  
  //build out lin list
  vector<vector<size_t>> equiv_linmom_combos=get_equiv_list(linmoms,"equiv_linmoms.txt");
  for(const auto &p : equiv_linmom_combos)
    out.linmoms.push_back({p[0]});
  
  //build out bil combo
  vector<vector<size_t>> equiv_bilmom_combos=get_equiv_list(bilmoms,"equiv_bilmoms.txt");
  fill_output_equivalent_momenta(out.bilmoms,equiv_linmom_combos,equiv_bilmom_combos,bilmoms);
  
  cout<<"Equiv bil:"<<endl;
  // for(auto &p : out.bilmoms)
  // for(auto &p : equiv_bilmom_combo)
  // {
  //   for(auto &pi : p)
  // 	cout<<pi<<endl;
  //   cout<<endl;
  // }
  
  out.set_indices();
  out.allocate();
  
  //average Zq
  for(size_t i=0;i<out.im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(i);
      const size_t out_ilinmom_combo=out_im_r_ilinmom_comp[2];
      
      for(const auto &t : get_Zq_tasks(out))
	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq : equiv_linmom_combos[out_ilinmom_combo]) ave+=(*t.in)[ieq];
  	  ave/=equiv_linmom_combos[out_ilinmom_combo].size();
  	}
    }
  
  //average pr_bil and Z
  for(size_t i=0;i<out.im_r_im_r_iZbil_ibilmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_iZbil_ibilmom_comp=out.im_r_im_r_iZbil_ibilmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_iZbil_ibilmom_comp[5];
      
      for(const auto &t : concat(get_pr_bil_tasks(out),get_Zbil_tasks(out)))
  	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq : equiv_bilmom_combos[out_imom_combo]) ave+=(*t.in)[ieq];
  	  ave/=equiv_bilmom_combos[out_imom_combo].size();
  	}
    }
  
  if(not Zbil_computed or recompute_Zbil) out.compute_Zbil();
  
  return out;
}

void ingredients_t::plot_Z(const string &suffix) const
{
  ingredients_t dum;
  
  //Zq
  for(auto &t : get_Zq_tasks(dum))
    {
      const djvec_t &Z=*t.in;
      const string &tag=t.tag;
      
      grace_file_t out("plots/"+tag+(suffix!=""?("_"+suffix):string(""))+".xmg");
      
      for(size_t im=0;im<_nm;im++)
  	for(size_t r=0;r<_nr;r++)
	  {
	    out.new_data_set();
	    for(size_t imom=0;imom<linmoms.size();imom++)
	      {
		const double p2hat=glb_moms[linmoms[imom][0]].p(L).tilde().norm2();
		out.write_ave_err(p2hat,Z[im_r_ilinmom_ind({im,r,imom})].ave_err());
	      }
	  }
    }
  
  //Zbil
  for(const auto &t : get_Zbil_tasks(dum))
    {
      //decript tuple
      const djvec_t &Z=*t.in;
      const string &tag=t.tag;
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  	{
	  grace_file_t out("plots/"+tag+"_"+Zbil_tag[iZbil]+(suffix!=""?("_"+suffix):string(""))+".xmg");
	  
	  //write mass by mass, only half of the combos
  	  for(size_t im1=0;im1<_nm;im1++)
	    for(size_t im2=im1;im2<_nm;im2++)
	      for(size_t r=0;r<_nr;r++)
		{
		  out.new_data_set();
		  
		  for(size_t imom=0;imom<linmoms.size();imom++)
		    {
		      const double p2hat=glb_moms[bilmoms[imom][0]].p(L).tilde().norm2();
		      out.write_ave_err(p2hat,Z[im_r_im_r_iZbil_ibilmom_ind({im1,r,im2,r,iZbil,imom})].ave_err());
		    }
		  out.new_data_set();
		}
	}
    }
}
