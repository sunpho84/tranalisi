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
  _nm=glb::nm;
  _nr=glb::nr;
  _am=glb::am;
  
  switch(scheme)
    {
    case RI_MOM:
      set_ri_mom_moms();
      break;
    }
}

void ingredients_t::set_ri_mom_moms()
{
  mom_combo.resize(glb_moms.size());
  for(size_t imom=0;imom<glb_moms.size();imom++) mom_combo[imom].push_back(imom);
}

void ingredients_t::ri_mom()
{
  vector<m_r_mom_conf_props_t> props; //!< store props for individual conf
  
  vector<task_list_t> read_tasks=prepare_read_prop_taks(props,conf_list);
  
  for(size_t imom=0;imom<mom_combo.size();imom++)
    {
      vector<jm_r_mom_props_t> jprops(glb::im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),use_QED); //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    size_t i_in_clust_hit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", hit "<<ihit+1<<"/"<<nhits<<", momentum "<<imom+1<<"/"<<mom_combo.size()<<endl;
	    read_time.start();
	    read_tasks[i_in_clust_hit].assolve_all(RECYCLE);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackknifed_props(jprops,props,use_QED,glb::im_r_ind,deltam_cr);
	    build_props_time.stop();
	    
	    build_verts_time.start();
	    build_all_mr_gbil_jackknifed_verts(jverts,props,im_r_im_r_igam_ind,im_r_ijack_ind,use_QED,deltam_cr);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackknifed_props(jprops,use_QED,clust_size);
      jverts.clusterize_all(use_QED,clust_size);
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
	  const size_t im_r_imom=im_r_imom_ind({im,r,imom});
	  
	  //compute inverse
	  invert_time.start();
	  prop_t prop_inv=jprops[im_r].LO[ijack].inverse();
	  //if(im_r_ijack==0) cout<<jprops[im_r].LO[ijack](0,0)<<endl;
	  jprop_inv[im_r][ijack]=prop_inv;
	  invert_time.stop();
	  
	  //compute Zq
	  Zq_time.start();
	  Zq[im_r_imom][ijack]=compute_Zq(prop_inv,imom);
	  Zq_sig1[im_r_imom][ijack]=compute_Zq_sig1(prop_inv,imom);
	  Zq_time.stop();
	  
	  //do the same with QED
	  if(use_QED)
	    {
	      invert_time.start();
	      prop_t prop_EM_inv=prop_inv*jprops[im_r].EM[ijack]*prop_inv;
	      jprop_EM_inv[im_r][ijack]=prop_EM_inv;
	      invert_time.stop();
	      
	      Zq_time.start();
	      Zq_sig1_EM[im_r_imom][ijack]=-compute_Zq_sig1(prop_EM_inv,imom);
	      Zq_time.stop();
	    }
	}
      
      //! an index running on all packed combo, and momenta
      const index_t all_imom_ind({{"All",im_r_im_r_iZbil_ind.max()},{"mom",mom_combo.size()}});
      
      proj_time.start();
      djvec_t pr_bil_temp=compute_proj_bil(jprop_inv,jverts.LO,jprop_inv,glb::im_r_ind);
      
      //QED
      djvec_t pr_bil_QED_temp;
      if(use_QED)
	{
	  const djvec_t pr_bil_EM=compute_proj_bil(jprop_inv,jverts.EM,jprop_inv,glb::im_r_ind);
	  const djvec_t pr_bil_a=compute_proj_bil(jprop_EM_inv,jverts.LO,jprop_inv,glb::im_r_ind);
	  const djvec_t pr_bil_b=compute_proj_bil(jprop_inv,jverts.LO,jprop_EM_inv,glb::im_r_ind);
	  pr_bil_QED_temp=pr_bil_a+pr_bil_b-pr_bil_EM;
	}
      
      //store
      vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil,&pr_bil_temp}};
      if(use_QED) store_task.push_back({&pr_bil_QED,&pr_bil_QED_temp});
      for(auto &t : store_task)
	for(size_t iall=0;iall<im_r_im_r_iZbil_ind.max();iall++)
	  (*t.first)[all_imom_ind({iall,imom})]=(*t.second)[iall];
      
      proj_time.stop();
    }
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
  imom_ind.set_ranges({{"mom",mom_combo.size()}});
  r_imom_ind.set_ranges({{"r",_nr},{"mom",mom_combo.size()}});
  im_r_imom_ind.set_ranges({{"m",_nm},{"r",_nr},{"mom",mom_combo.size()}});
  // indep_imom_ind.set_ranges({{"indep_mom",equiv_imoms.size()}});
  // r_indep_imom_ind.set_ranges({{"r",nr},{"indep_mom",equiv_imoms.size()}});
  // im_r_indep_imom_ind.set_ranges({{"m",nm},{"r",nr},{"indep_mom",equiv_imoms.size()}});
  iZbil_imom_ind.set_ranges({{"iZbil",nZbil},{"mom",mom_combo.size()}});
  im_r_im_r_iZbil_imom_ind=im_r_im_r_iZbil_ind*index_t({{"mom",mom_combo.size()}});
}

void ingredients_t::allocate()
{
  Zq.resize(im_r_imom_ind.max());
  Zq_sig1.resize(im_r_imom_ind.max());
  if(use_QED) Zq_sig1_EM.resize(im_r_imom_ind.max());
  
  pr_bil.resize(im_r_im_r_iZbil_imom_ind.max());
  if(use_QED) pr_bil_QED.resize(im_r_im_r_iZbil_imom_ind.max());
  
  Zbil.resize(im_r_im_r_iZbil_imom_ind.max());
  if(use_QED) Zbil_QED.resize(im_r_im_r_iZbil_imom_ind.max());
}

void ingredients_t::create_from_scratch(const string ingredients_path)
{
  set_pars_for_scratch();
  
  set_indices();
  
  allocate();
  
  if(file_exists(ingredients_path)) bin_read(ingredients_path);
  else
    {
      switch(scheme)
	{
	case RI_MOM:
	  ri_mom();
	  break;
	}
      bin_write(ingredients_path);
    }
  
  compute_Zbil();
}

void ingredients_t::compute_Zbil()
{
  Zbil_computed=true;
  
  for(size_t imom=0;imom<mom_combo.size();imom++)
    for(size_t im_r_im_r_iZbil=0;im_r_im_r_iZbil<im_r_im_r_iZbil_ind.max();im_r_im_r_iZbil++)
      {
	const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
	const vector<size_t> im_r1_comp=subset(im_r_im_r_iZbil_comp,0,2);
	const vector<size_t> im_r2_comp=subset(im_r_im_r_iZbil_comp,2,4);
	const size_t iZbil=im_r_im_r_iZbil_comp[4];
	const size_t im_r1_imom=im_r_imom_ind(concat(im_r1_comp,imom));
	const size_t im_r2_imom=im_r_imom_ind(concat(im_r2_comp,imom));
	
	const size_t im_r_im_r_iZbil_imom=im_r_im_r_iZbil_imom_ind(concat(im_r1_comp,im_r2_comp,vector<size_t>({iZbil,imom})));
	
	Zbil[im_r_im_r_iZbil_imom]=
	  sqrt(Zq_sig1[im_r1_imom]*Zq_sig1[im_r2_imom])/pr_bil[im_r_im_r_iZbil_imom];
	
	if(use_QED)
	  {
	    Zbil_QED[im_r_im_r_iZbil_imom]=
	      pr_bil_QED[im_r_im_r_iZbil_imom]/pr_bil[im_r_im_r_iZbil_imom]+
	      (Zq_sig1_EM[im_r1_imom]/Zq_sig1[im_r1_imom]+Zq_sig1_EM[im_r2_imom]/
	       Zq_sig1[im_r2_imom])/2.0;
	  }
      }
}

ingredients_t ingredients_t::chir_extrap() const
{
  ingredients_t out;
  out._nr=1;
  out._nm=1;
  out._am={0.0};
  out.mom_combo=mom_combo;
  out.Zbil_computed=Zbil_computed;
  
  out.set_indices();
  out.allocate();
  
  //! list of task to chirally extrapolate Zq
  vector<tuple<const djvec_t*,djvec_t*,string>> Zq_chirextr_tasks{
    {&Zq,&out.Zq,string("Zq")},
    {&Zq_sig1,&out.Zq_sig1,"Zq_sig1"}};
  if(use_QED) Zq_chirextr_tasks.push_back(make_tuple(&Zq_sig1_EM,&out.Zq_sig1_EM,"Zq_sig1_EM"));
  
  //! list of tasks to chirally extrapolate bilinears
  vector<tuple<const djvec_t*,djvec_t*,string>>
    pr_bil_chirextr_tasks{{&pr_bil,&out.pr_bil,string("pr_bil")}};
  if(use_QED) pr_bil_chirextr_tasks.push_back(make_tuple(&pr_bil_QED,&out.pr_bil_QED,string("pr_bil_QED")));
  
  for(size_t imom_combo=0;imom_combo<mom_combo.size();imom_combo++)
    {
      const size_t imom=mom_combo[imom_combo][0];
      
      //extrapolate to chiral limit Zq
      for(auto & p : Zq_chirextr_tasks)
  	{
  	  djvec_t y(_nm);
  	  const djvec_t &Zq=(*get<0>(p));
  	  djvec_t &Zq_chir=(*get<1>(p));
  	  const string &tag=get<2>(p);
	  
  	  //slice m
  	  double am_max=*max_element(_am.begin(),_am.end())*1.1;
  	  for(size_t im=0;im<_nm;im++)
  	    {
  	      //averages r if both asked
  	      y[im]=0.0;
  	      for(size_t r=0;r<_nr;r++)
  		{
  		  size_t i=im_r_imom_ind({im,r,imom_combo});
  		  y[im]+=Zq[i]/_nr;
  		  cout<<tag<<"["<<i<<"=(im"<<im<<"r"<<r<<"imom_combo"<<imom_combo<<")], mom "<<glb_moms[imom].p(L).norm2()<<": "<<Zq[i].ave_err()<<endl;
  		}
  	      cout<<tag<<"[im"<<im<<"imom"<<imom<<"]: "<<y[im][0]<<" "<<y[im].err()<<endl;
  	    }
	  
  	  //fit and write the result
  	  djvec_t coeffs=poly_fit(_am,y,1,am_min,am_max);
  	  if(imom%print_each_mom==0)
  	    {
  	      grace_file_t plot("plots/chir_extr_"+tag+"_mom_"+to_string(imom)+".xmg");
  	      write_fit_plot(plot,0,am_max,bind(poly_eval<djvec_t>,coeffs,_1),_am,y);
  	      plot.write_ave_err(0,coeffs[0].ave_err());
  	    }
  	  //extrapolated value
  	  Zq_chir[imom]=coeffs[0];
  	}
      
      //extrapolate to chiral limit bil
      for(auto & p : pr_bil_chirextr_tasks)
	for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	  {
	    const djvec_t &pr=(*get<0>(p));
	    djvec_t &pr_chir=(*get<1>(p));
	    const string &tag=get<2>(p);
	    
	    //check if we need to subtract the pole
	    const bool sub_pole=(iZbil==iZS or iZbil==iZP);
	    
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
  		  y_plot[i]=0.0;
  		  for(size_t r=0;r<_nr;r++) y_plot[i]+=pr[im_r_im_r_iZbil_imom_ind({im1,r,im2,r,iZbil,imom})]/_nr;
		  
  		  if(sub_pole) y[i]=x[i]*y_plot[i];
  		  else         y[i]=y_plot[i];
  		  //increment
  		  i++;
  		}
	    
  	    //fit and store extrapolated value
	    const size_t iout=out.im_r_im_r_iZbil_imom_ind({0,0,0,0,iZbil,imom});
  	    const djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1),2.0*am_min,2.0*am_max);
  	    pr_chir[iout]=coeffs[sub_pole?1:0];
	    
  	    //plot
  	    if(imom%print_each_mom==0)
  	      {
  		grace_file_t plot("plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_mom_"+to_string(imom)+".xmg");
  		write_fit_plot(plot,2*am_min,2*am_max,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
  		plot.write_ave_err(0.0,pr_chir[iout].ave_err());
  	      }
  	  }
    }
  
  /////////////////////////////////////////////
  
  out.compute_Zbil();
  
  return out;
}

ingredients_t ingredients_t::subtract_Oa2(const bool recompute_Zbil) const
{
  ingredients_t out=(*this);
  
  //subtract Zq
  for(size_t i=0;i<im_r_imom_ind.max();i++)
    {
      const vector<size_t> im_r_imom_comp=im_r_imom_ind(i);
      const size_t imom_combo=im_r_imom_comp[2];
      const size_t imom=mom_combo[imom_combo][0];
      
      imom_t mom=glb_moms[imom];
      double sub=g2tilde*sig1_a2(act,gf::LANDAU,group::SU3,mom,L);
      out.Zq[i]=Zq[i]-sub;
      out.Zq_sig1[i]=Zq_sig1[imom]-sub;
      if(use_QED)
	{
	  double sub_EM=-1.0*/*coupling gets -1 due to definition of expansion*/sig1_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L);
	  out.Zq_sig1_EM[i]=Zq_sig1_EM[i]-sub_EM*Zq_sig1[i]; //factorization hypotesis
	}
    }
  
  //subtract bil
  for(size_t i=0;i<im_r_im_r_iZbil_imom_ind.max();i++)
    {
      const vector<size_t> im_r_im_r_iZbil_imom_comp=im_r_im_r_iZbil_imom_ind(i);
      const size_t iZbil=im_r_im_r_iZbil_imom_comp[4];
      const size_t imom_combo=im_r_im_r_iZbil_imom_comp[5];
      const size_t imom=mom_combo[imom_combo][0];
      
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
  for(size_t i=0;i<im_r_imom_ind.max();i++)
    {
      const vector<size_t> im_r_imom_comp=im_r_imom_ind(i);
      const size_t imom_combo=im_r_imom_comp[2];
      const size_t imom=mom_combo[imom_combo][0];
      
      imom_t mom=glb_moms[imom];
      double p2=mom.p(L).norm2();
      double evolver_Zq=evolution_Zq_to_RIp(Nf,ord,ainv,p2);
      cout<<"EvolverZq["<<p2<<"]="<<evolver_Zq<<endl;
      
      out.Zq[i]=Zq[i]/evolver_Zq;
      out.Zq_sig1[i]=Zq_sig1[imom]/evolver_Zq;
      if(use_QED) out.Zq_sig1_EM[i]=Zq_sig1_EM[i]/evolver_Zq; //neglecting QED evolution
    }
  
  if(not Zbil_computed) CRASH("Cannot evolve before computing Zbil");
  
  //evolve bil
  for(size_t i=0;i<im_r_im_r_iZbil_imom_ind.max();i++)
    {
      const vector<size_t> im_r_im_r_iZbil_imom_comp=im_r_im_r_iZbil_imom_ind(i);
      const size_t iZbil=im_r_im_r_iZbil_imom_comp[4];
      const size_t imom_combo=im_r_im_r_iZbil_imom_comp[5];
      const size_t imom=mom_combo[imom_combo][0];
      
      const imom_t mom=glb_moms[imom];
      const double p2=mom.p(L).norm2();
      const double evolver_Zbil=evolution_Zbil_to_RIp(iZbil_t_list[iZbil],Nf,ord,ainv,p2);
      cout<<"EvolverZ"<<Zbil_tag[iZbil]<<"["<<p2<<"]="<<evolver_Zbil<<endl;
      
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
  
  const size_t np=mom_combo[0].size(); //number of momenta defining the combo
  
  //periodicity
  enum pe_t{UNK,PER,APE};
  const string pe_tag[3]={"UNK","PE","APE"};
  pe_t pe=UNK;
  double p0=fabs(ph_mom[0]);
  if(fabs(p0)<1e-10) pe=PER;
  if(fabs(p0-0.5)<1e-10) pe=APE;
  cout<<"Periodicity: "<<pe_tag[pe]<<endl;
  
  map<vector<imom_t>,vector<size_t>> equiv_mom_combo_map;
  for(size_t imom_combo=0;imom_combo<mom_combo.size();imom_combo++)
    {
      //get representative
      vector<imom_t> cr(np);
      for(size_t ip=0;ip<np;ip++)
	{
	  const size_t imom=mom_combo[imom_combo][ip];
	  //decide time component
	  cr[ip][0]=glb_moms[imom][0];
	  if(cr[ip][0]<0)
	    switch(pe)
	      {
	      case UNK: CRASH("phase on momentum 0 cannot be %lg",ph_mom[0]);break;
	      case PER: cr[ip][0]=-cr[ip][0];break;
	      case APE: cr[ip][0]=-cr[ip][0]-1;break;
	      }
	  
	  //decide space components
	  for(size_t mu=1;mu<NDIM;mu++) cr[ip][mu]=abs(glb_moms[imom][mu]);
	  sort(&cr[ip][1],cr[ip].end());
	}
      sort(cr.begin(),cr.end());
      //store the index to equvalents
      equiv_mom_combo_map[cr].push_back(imom_combo);
    }
  
  //prepare the index of the output
  const size_t nequiv=equiv_mom_combo_map.size();
  out.mom_combo.reserve(nequiv);
  
  //store the output vector to which the input contributes
  vector<vector<size_t>> equiv_of_combo;
  
  //fill the output index and write it to file
  ofstream mom_out("equiv_mom_combo.txt");
  mom_out<<"Found "<<nequiv<<" independent momenta combo"<<endl;
  for(auto &combo_class : equiv_mom_combo_map)
    {
      //add to the output list of momenta
      const size_t combo_repr=combo_class.second[0]; //take the first real combo as the key could not exist
      out.mom_combo.push_back(mom_combo[combo_repr]);
      const size_t imom=mom_combo[combo_repr][0];
      
      //trasform map into vector
      equiv_of_combo.push_back(combo_class.second);
      
      //print details
      mom_out<<imom<<" , p2hat: "<<glb_moms[imom].p(L).tilde().norm2()<<endl;
      mom_out<<" Equivalent to: "<<combo_class.second.size()<<" mom combos: "<<endl;
      for(auto &eq_combo : combo_class.second)
  	{
	  for(size_t ip=0;ip<np;ip++)
	    {
	      const size_t imom=mom_combo[eq_combo][ip];
	      mom_out<<"  "<<imom<<"=";
	      //components
	      mom_out<<"{";
	      for(size_t mu=0;mu<NDIM;mu++)
		{
		  if(mu) mom_out<<",";
		  mom_out<<glb_moms[imom][mu];
		}
	      mom_out<<"}";
	    }
	  mom_out<<endl;
	}
    }
  
  out.set_indices();
  out.allocate();
  
  //average Zq
  vector<tuple<const djvec_t*,djvec_t*>>
    Zq_ave_tasks{{&Zq,&out.Zq},{&Zq_sig1,&out.Zq_sig1}};
  if(use_QED) Zq_ave_tasks.push_back({&Zq_sig1_EM,&out.Zq_sig1_EM});
  for(size_t i=0;i<out.im_r_imom_ind.max();i++)
    {
      const vector<size_t> out_im_r_imom_comp=out.im_r_imom_ind(i);
      const size_t out_imom_combo=out_im_r_imom_comp[2];
      
      for(const auto &p : Zq_ave_tasks)
	{
	  djack_t &ave=(*get<1>(p))[i];
	  ave=0.0;
	  for(const size_t ieq : equiv_of_combo[out_imom_combo]) ave+=(*get<0>(p))[ieq];
	  ave/=equiv_of_combo[out_imom_combo].size();
	}
    }
  
  //average pr_bil
  vector<tuple<const djvec_t*,djvec_t*>>
    bil_ave_tasks{{&pr_bil,&out.pr_bil},{&Zbil,&out.Zbil}};
  if(use_QED)
    {
      bil_ave_tasks.push_back({&pr_bil_QED,&out.pr_bil_QED});
      bil_ave_tasks.push_back({&Zbil_QED,&out.Zbil_QED});
    }
  for(size_t i=0;i<out.im_r_im_r_iZbil_imom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_iZbil_imom_comp=out.im_r_im_r_iZbil_imom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_iZbil_imom_comp[5];
      
      for(const auto &p : bil_ave_tasks)
	{
	  djack_t &ave=(*get<1>(p))[i];
	  ave=0.0;
	  for(const size_t ieq : equiv_of_combo[out_imom_combo]) ave+=(*get<0>(p))[ieq];
	  ave/=equiv_of_combo[out_imom_combo].size();
	}
    }
  
  if(not Zbil_computed or recompute_Zbil) out.compute_Zbil();
  
  return out;
}

void ingredients_t::plot_Z(const string &suffix) const
{
  vector<tuple<const djvec_t*,string>> Zq_tasks={{&Zq,"Zq"},{&Zq_sig1,"Zq_sig1"}};
  if(use_QED) Zq_tasks.push_back(make_tuple(&Zq_sig1_EM,"Zq_sig1_EM"));
  for(const auto &p : Zq_tasks)
    {
      const djvec_t &Z=*get<0>(p);
      const string &tag=get<1>(p);
      grace_file_t out("plots/"+tag+(suffix!=""?("_"+suffix):string(""))+".xmg");
      
      for(size_t i=0;i<im_r_imom_ind.max();i++)
	{
	  const vector<size_t> im_r_imom_comp=im_r_imom_ind(i);
	  const size_t imom_combo=im_r_imom_comp[2];
	  const size_t imom=mom_combo[imom_combo][0];
	  
	  const double p2hat=glb_moms[imom].p(L).tilde().norm2();
	  out.write_ave_err(p2hat,Z[i].ave_err());
	}
    }
}
