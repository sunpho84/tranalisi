#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <MOM2/Zbil.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

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

djvec_t perens_t::compute_proj_bil(const vjqprop_t &jprop_inv_in,const vector<jqprop_t> &jverts,const vjqprop_t &jprop_inv_ou) const
{
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
	  const size_t im_r_im_r_iG=im_r_im_r_igam_ind(im_r_im_r_iG_comp);
	  
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

void perens_t::build_all_mr_gbil_jackkniffed_verts(jbil_vert_t &jbil,const vector<m_r_mom_conf_qprops_t> &props_in,const vector<m_r_mom_conf_qprops_t> &props_ou) const
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
      const m_r_mom_conf_qprops_t &p_in=props_in[im_r_ijack_ind({im_fw,r_fw,iclust})];
      const m_r_mom_conf_qprops_t &p_ou=props_ou[im_r_ijack_ind({im_bw,r_bw,iclust})];
      
      //create list of operations
      vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>> list={{&jbil.LO[im_r_im_r_igam],&p_in.LO,&p_ou.LO}};
      if(pars::use_QED)
	for(auto &o : vector<tuple<jqprop_t*,const qprop_t*,const qprop_t*>>({
	      {&jbil.PH[im_r_im_r_igam],&p_in.F,&p_ou.F},
	      {&jbil.PH[im_r_im_r_igam],&p_in.FF,&p_ou.LO},
	      {&jbil.PH[im_r_im_r_igam],&p_in.LO,&p_ou.FF},
	      {&jbil.PH[im_r_im_r_igam],&p_in.T,&p_ou.LO},
	      {&jbil.PH[im_r_im_r_igam],&p_in.LO,&p_ou.T},
	      {&jbil.CR_CT_in[im_r_im_r_igam],&p_in.P,&p_ou.LO},
	      {&jbil.CR_CT_ou[im_r_im_r_igam],&p_in.LO,&p_ou.P},
	      {&jbil.TM_CT_in[im_r_im_r_igam],&p_in.S,&p_ou.LO},
	      {&jbil.TM_CT_ou[im_r_im_r_igam],&p_in.LO,&p_ou.S}
	    }))
	  list.push_back(o);
      
      //create the vertex
      for(auto &o : list)
	build_jackkniffed_vert_Gamma(*get<0>(o),*get<1>(o),iG,*get<2>(o),iclust);
    }
}

void perens_t::mom_compute_bil()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    {
      const size_t imom1=bilmoms[ibilmom][1];
      const size_t imom2=bilmoms[ibilmom][2];
      const bool read2=(imom1!=imom2);
      
      vector<jm_r_mom_qprops_t> jprops1(im_r_ind.max()); //!< jackknived props
      vector<jm_r_mom_qprops_t> jprops2(im_r_ind.max()); //!< jackknived props
      jbil_vert_t jverts(im_r_im_r_igam_ind.max(),pars::use_QED);   //!< jackknived vertex
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t mom1=linmoms[imom1][0];
	    const size_t mom2=linmoms[imom2][0];
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    cout<<"Working on qbil, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum combo "<<ibilmom+1<<"/"<<bilmoms.size()<<", "
	      "moms: "<<mom1<<" "<<mom2<<endl;
	    
	    //read
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props1=read_all_qprops_mom(files,i_in_clust_ihit,mom1);
	    const vector<m_r_mom_conf_qprops_t> props2=(read2?read_all_qprops_mom(files,i_in_clust_ihit,mom2):props1);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops1,props1);
	    build_all_mr_jackkniffed_qprops(jprops2,props2);
	    build_props_time.stop();
	    
	    //build all bilinear verts
	    build_verts_time.start();
	    build_all_mr_gbil_jackkniffed_verts(jverts,props1,props2);
	    build_verts_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops1);
      clusterize_all_mr_jackkniffed_qprops(jprops2);
      jverts.clusterize_all(pars::use_QED,clust_size,im_r_im_r_igam_ind,deltam_cr,deltam_tm);
      clust_time.stop();
      
      vector<jqprop_t> jprop_inv1; //!< inverse propagator1
      vector<jqprop_t> jprop_inv2; //!< inverse propagator2
      vector<jqprop_t> jprop_QED_inv1; //!< inverse propagator1 with em insertion
      vector<jqprop_t> jprop_QED_inv2; //!< inverse propagator2 with em insertion
      
      get_inverse_propagators(jprop_inv1,jprop_QED_inv1,jprops1);
      get_inverse_propagators(jprop_inv2,jprop_QED_inv2,jprops2);
      
      proj_time.start();
      djvec_t pr_bil_temp=compute_proj_bil(jprop_inv1,jverts.LO,jprop_inv2);
      
      //QED
      djvec_t pr_bil_QED_temp;
      if(pars::use_QED)
	{
	  const djvec_t pr_bil_QED=compute_proj_bil(jprop_inv1,jverts.QED,jprop_inv2);
	  const djvec_t pr_bil_a=compute_proj_bil(jprop_QED_inv1,jverts.LO,jprop_inv2);
	  const djvec_t pr_bil_b=compute_proj_bil(jprop_inv1,jverts.LO,jprop_QED_inv2);
	  pr_bil_QED_temp=
	    (-pr_bil_a
	     -pr_bil_b
	     +pr_bil_QED);
	}
      
      //! an index running on all packed combo, and momenta
      const index_t all_ibilmom_ind({{"All",im_r_im_r_iZbil_ind.max()},{"bilmom",bilmoms.size()}});
      
      //store
      vector<pair<djvec_t*,djvec_t*>> store_task={{&pr_bil,&pr_bil_temp}};
      if(pars::use_QED) store_task.push_back({&pr_bil_QED,&pr_bil_QED_temp});
      for(auto &t : store_task)
	for(size_t iall=0;iall<im_r_im_r_iZbil_ind.max();iall++)
	  (*t.first)[all_ibilmom_ind({iall,ibilmom})]=(*t.second)[iall];
      
      proj_time.stop();
    }
}

vector<perens_t::task_t> perens_t::get_pr_bil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_pr_bil,in_pr_bil_QED;
  for(auto &e : ens)
    {
      in_pr_bil.push_back(&e->pr_bil);
      if(pars::use_QED) in_pr_bil_QED.push_back(&e->pr_bil_QED);
    }
  vector<task_t> pr_bil_tasks={{&pr_bil,in_pr_bil,im_r_im_r_iZbil_ibilmom_ind,"pr_bil",QCD_task}};
  if(pars::use_QED) pr_bil_tasks.push_back({&pr_bil_QED,in_pr_bil_QED,im_r_im_r_iZbil_ibilmom_ind,"pr_bil_QED",QED_task});
  
  return pr_bil_tasks;
}

vector<perens_t::task_t> perens_t::get_Zbil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_Zbil,in_Zbil_QED_rel;
  for(auto &e : ens)
    {
      in_Zbil.push_back(&e->Zbil);
      if(pars::use_QED) in_Zbil_QED_rel.push_back(&e->Zbil_QED_rel);
    }
  vector<task_t> Zbil_tasks={{&Zbil,in_Zbil,im_r_im_r_iZbil_ibilmom_ind,"Zbil",QCD_task}};
  if(pars::use_QED) Zbil_tasks.push_back({&Zbil_QED_rel,in_Zbil_QED_rel,im_r_im_r_iZbil_ibilmom_ind,"Zbil_QED_rel",QED_task});
  
  return Zbil_tasks;
}

void perens_t::compute_Zbil()
{
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
	
	if(pars::use_QED)
	    Zbil_QED_rel[im_r_im_r_iZbil_ibilmom]=
	      -pr_bil_QED[im_r_im_r_iZbil_ibilmom]/pr_bil[im_r_im_r_iZbil_ibilmom]
	      +(Zq_sig1_QED[im_r1_ilinmom1]/Zq_sig1[im_r1_ilinmom1]+
		Zq_sig1_QED[im_r2_ilinmom2]/Zq_sig1[im_r2_ilinmom2])/2.0;
      }
}

void perens_t::plot_Zbil(const string &suffix)
{
  for(const auto &t : this->get_Zbil_tasks())
    {
      //decript tuple
      const djvec_t &Z=*t.out;
      const string &tag=t.tag;
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  	{
	  grace_file_t out(dir_path+"/plots/"+tag+"_"+Zbil_tag[iZbil]+(suffix!=""?("_"+suffix):string(""))+".xmg");
	  
	  //write mass by mass, only half of the combos
  	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=im1;im2<nm;im2++)
	      for(size_t r=0;r<nr;r++)
		{
		  out.new_data_set();
		  
		  for(size_t imom=0;imom<bilmoms.size();imom++)
		    {
		      const double p2tilde=all_moms[bilmoms[imom][0]].p(L).tilde().norm2();
		      out.write_ave_err(p2tilde,Z[im_r_im_r_iZbil_ibilmom_ind({im1,r,im2,r,iZbil,imom})].ave_err());
		    }
		}
	}
    }
}

void perens_t::average_r_Zbil(perens_t &out) const
{
  cout<<"Averaging r for bil"<<endl;
  
  for(auto &t : out.get_bil_tasks({this}))
    {
      cout<<" "<<t.tag<<endl;
      
      const djvec_t &pr=*t.in.front();
      djvec_t &pr_rave=*t.out;
      
      for(size_t out_i=0;out_i<out.im_r_im_r_iZbil_ibilmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_im_r_iZbil_ibilmom_comp=out.im_r_im_r_iZbil_ibilmom_ind(out_i);
	  vector<size_t> im_r_im_r_iZbil_ibilmom_comp=out_im_r_im_r_iZbil_ibilmom_comp;
	  
	  pr_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_im_r_iZbil_ibilmom_comp[1]=
		im_r_im_r_iZbil_ibilmom_comp[3]=r;
	      const size_t i=im_r_im_r_iZbil_ibilmom_ind(im_r_im_r_iZbil_ibilmom_comp);
	      pr_rave[out_i]+=pr[i];
	    }
	  pr_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_Zbil(perens_t &out,const vector<vector<size_t>> &equiv_bilmom_combos) const
{
  for(size_t i=0;i<out.im_r_im_r_iZbil_ibilmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_im_r_iZbil_ibilmom_comp=out.im_r_im_r_iZbil_ibilmom_ind(i);
      const size_t out_imom_combo=out_im_r_im_r_iZbil_ibilmom_comp[5];
      
      for(const auto &t : out.get_bil_tasks({this}))
  	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq_mom : equiv_bilmom_combos[out_imom_combo])
	    {
	      vector<size_t> im_r_im_r_iZbil_ibilmom_comp=out_im_r_im_r_iZbil_ibilmom_comp;
	      im_r_im_r_iZbil_ibilmom_comp[5]=ieq_mom;
	      const size_t ieq=im_r_im_r_iZbil_ibilmom_ind(im_r_im_r_iZbil_ibilmom_comp);
	      ave+=(*t.in.front())[ieq];
	    }
  	  ave/=equiv_bilmom_combos[out_imom_combo].size();
  	}
    }
}

void perens_t::val_chir_extrap_Zbil(perens_t &out) const
{
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(auto &t : out.get_pr_bil_tasks({this}))
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  const djvec_t &pr=*t.in.front();
	  djvec_t &pr_chir=*t.out;
	  const string &tag=t.tag;
	  
	  //check if we need to subtract the pole
	  const bool sub_pole=(iZbil==iZS or iZbil==iZP);
	  const size_t x_pow=(sub_pole?
			      (t.QCD_QED_task==QED_task?2:1)
			      :0);
	  
	  //open the plot file if needed
	  const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_bilmom_"+to_string(ibilmom)+".xmg";
	  grace_file_t *plot=nullptr;
	  if(ibilmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	  
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
		      
		      if(std::isnan(x[i])) CRASH("Nanning %d",i);
		      
		      //compute y and y_plot
		      y_plot[i]=pr[im_r_im_r_iZbil_ibilmom_ind({im1,r1,im2,r2,iZbil,ibilmom})];
		      //fit x*y if pole present
		      y[i]=pow(x[i],x_pow)*y_plot[i];
		      //increment the number of mass combos
		      i++;
		    }
		
		//fit, store and write the result
		const djvec_t coeffs=poly_fit(x,y,1+x_pow);
		const size_t iout=out.im_r_im_r_iZbil_ibilmom_ind({0,r1,0,r2,iZbil,ibilmom});
		pr_chir[iout]=coeffs[x_pow];
		if(plot!=nullptr)
		  {
		    auto xminmax=minmax_element(x.begin(),x.end());
		    double xmin=*xminmax.first*0.5;
		    double xmax=*xminmax.second*1.1;
		    write_fit_plot(*plot,xmin,xmax,[&coeffs,sub_pole,x_pow](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/pow(x,x_pow);},x,y_plot);
		    plot->write_ave_err(0.0,pr_chir[iout].ave_err());
		  }
	      }
	  if(plot) delete plot;
	}
}
