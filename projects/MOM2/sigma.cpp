#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <perens.hpp>
#include <timings.hpp>

void perens_t::plot_sigma(const string &suffix)
{
  cout<<"Plotting all sigma of "<<dir_path<<" for suffix: \""<<suffix<<"\""<<endl;
  
  auto tasks=this->get_sigma_tasks();
  
  for(auto &t : tasks)
    {
      const djvec_t &Z=*t.out;
      const string &tag=t.tag;
      
      grace_file_t out(dir_path+"/plots/"+tag+(suffix!=""?("_"+suffix):string(""))+".xmg");
      
      for(size_t im=0;im<nm;im++)
  	for(size_t r=0;r<nr;r++)
	  {
	    out.new_data_set();
	    for(size_t imom=0;imom<linmoms.size();imom++)
	      {
		const double p2tilde=all_moms[linmoms[imom][0]].p(L).tilde().norm2();
		out.write_ave_err(p2tilde,Z[im_r_ilinmom_ind({im,r,imom})].ave_err());
	      }
	  }
    }
}

vector<perens_t::task_t> perens_t::get_sigma_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_sigma1_LO,in_sigma2_LO,in_sigma3_LO;
  vector<const djvec_t*> in_sigma1_CR_CT,in_sigma2_CR_CT,in_sigma3_CR_CT;
  vector<const djvec_t*> in_sigma1_TM_CT,in_sigma2_TM_CT,in_sigma3_TM_CT;
  vector<const djvec_t*> in_sigma1_PH,in_sigma2_PH,in_sigma3_PH;
  
  for(auto &e : ens)
    {
      in_sigma1_LO.push_back(&e->sigma1_LO);
      in_sigma2_LO.push_back(&e->sigma2_LO);
      in_sigma3_LO.push_back(&e->sigma3_LO);
      ///
      if(pars::use_QED)
	{
	  in_sigma1_CR_CT.push_back(&e->sigma1_CR_CT);
	  in_sigma2_CR_CT.push_back(&e->sigma2_CR_CT);
	  in_sigma3_CR_CT.push_back(&e->sigma3_CR_CT);
	  ///
	  in_sigma1_TM_CT.push_back(&e->sigma1_TM_CT);
	  in_sigma2_TM_CT.push_back(&e->sigma2_TM_CT);
	  in_sigma3_TM_CT.push_back(&e->sigma3_TM_CT);
	  ///
	  in_sigma1_PH.push_back(&e->sigma1_PH);
	  in_sigma2_PH.push_back(&e->sigma2_PH);
	  in_sigma3_PH.push_back(&e->sigma3_PH);
	}
    }
  
  vector<task_t> sigma_tasks={
    {&sigma1_LO,in_sigma1_LO,im_r_ilinmom_ind,"sigma1_LO",QCD_task},
    {&sigma2_LO,in_sigma2_LO,im_r_ilinmom_ind,"sigma2_LO",QCD_task},
    {&sigma3_LO,in_sigma3_LO,im_r_ilinmom_ind,"sigma3_LO",QCD_task}};
  ///
  if(pars::use_QED)
    {
      sigma_tasks.push_back({&sigma1_CR_CT,in_sigma1_CR_CT,im_r_ilinmom_ind,"sigma1_CR_CT",QED_task});
      sigma_tasks.push_back({&sigma2_CR_CT,in_sigma2_CR_CT,im_r_ilinmom_ind,"sigma2_CR_CT",QED_task});
      sigma_tasks.push_back({&sigma3_CR_CT,in_sigma3_CR_CT,im_r_ilinmom_ind,"sigma3_CR_CT",QED_task});
      ///
      sigma_tasks.push_back({&sigma1_TM_CT,in_sigma1_TM_CT,im_r_ilinmom_ind,"sigma1_TM_CT",QED_task});
      sigma_tasks.push_back({&sigma2_TM_CT,in_sigma2_TM_CT,im_r_ilinmom_ind,"sigma2_TM_CT",QED_task});
      sigma_tasks.push_back({&sigma3_TM_CT,in_sigma3_TM_CT,im_r_ilinmom_ind,"sigma3_TM_CT",QED_task});
      ///
      sigma_tasks.push_back({&sigma1_PH,in_sigma1_PH,im_r_ilinmom_ind,"sigma1_PH",QED_task});
      sigma_tasks.push_back({&sigma2_PH,in_sigma2_PH,im_r_ilinmom_ind,"sigma2_PH",QED_task});
      sigma_tasks.push_back({&sigma3_PH,in_sigma3_PH,im_r_ilinmom_ind,"sigma3_PH",QED_task});
    }
  
  return sigma_tasks;
}

perens_t& perens_t::compute_sigmas()
{
  vector<raw_file_t> files=setup_read_all_qprops_mom(conf_list);
  
  for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
    {
      const size_t mom=linmoms[ilinmom][0];
      vector<jm_r_mom_qprops_t> jprops(im_r_ind.max()); //!< jackknived props
      
      for(size_t i_in_clust=0;i_in_clust<clust_size;i_in_clust++)
	for(size_t ihit=0;ihit<nhits_to_use;ihit++)
	  {
	    const size_t i_in_clust_ihit=i_in_clust_ihit_ind({i_in_clust,ihit});
	    const size_t mom=linmoms[ilinmom][0];
	    cout<<"Working on qprop, "
	      "clust_entry "<<i_in_clust+1<<"/"<<clust_size<<", "
	      "hit "<<ihit+1<<"/"<<nhits<<", "
	      "momentum "<<ilinmom+1<<"/"<<linmoms.size()<<", "
	      "mom: "<<mom<<endl;
	    read_time.start();
	    const vector<m_r_mom_conf_qprops_t> props=read_all_qprops_mom(files,i_in_clust_ihit,mom);
	    read_time.stop();
	    
	    //build all props
	    build_props_time.start();
	    build_all_mr_jackkniffed_qprops(jprops,props);
	    build_props_time.stop();
	  }
      
      //clusterize
      clust_time.start();
      clusterize_all_mr_jackkniffed_qprops(jprops);
      clust_time.stop();
      
      vector<jm_r_mom_qprops_t> jprops_inv(im_r_ind.max()); //!< jackknived inverse props
      
      get_inverse_propagators(jprops_inv,jprops);
      
#pragma omp parallel for reduction(+:invert_time,sigma_time)
      for(size_t im_r_ijack=0;im_r_ijack<im_r_ijackp1_ind.max();im_r_ijack++)
	{
	  //decript indices
	  const vector<size_t> im_r_ijack_comps=im_r_ijackp1_ind(im_r_ijack);
	  const size_t im=im_r_ijack_comps[0],r=im_r_ijack_comps[1],ijack=im_r_ijack_comps[2];
	  const size_t im_r=im_r_ind({im,r});
	  const size_t im_r_ilinmom=im_r_ilinmom_ind({im,r,ilinmom});
	  
	  const p_t ptilde=all_moms[mom].p(L).tilde();
	  
	  //! function to compute sigma1
	  auto compute_sigma1=[&](const jqprop_t prop_inv)
	    {
	      double out=0.0;
	      
	      for(size_t mu=0;mu<NDIM;mu++)
		if(fabs(ptilde[mu])>1e-10)
		  out+=
		    (prop_inv[ijack]*quaGamma[igmu[mu]]).trace().imag()/
		    (12.0*ptilde[mu]*V*all_moms[mom].Np());
	      
	      return out;
	    };
	  
	  //! function to compute sigma2
	  auto compute_sigma2=[&](const jqprop_t prop_inv)
	    {
	      return (prop_inv[ijack]*quaGamma[0]).trace().real()/(12.0*V);
	    };
	  
	  //! function to compute sigma3
	  auto compute_sigma3=[&](const jqprop_t prop_inv)
	    {
	      return (prop_inv[ijack]*quaGamma[5]).trace().imag()/(12.0*V);
	    };
	  
	  //! define the computation of the sigmas
#define COMPUTE_SIGMA(A)						\
	  sigma1_ ## A[im_r_ilinmom][ijack]=compute_sigma1(j.A);	\
	  sigma2_ ## A[im_r_ilinmom][ijack]=compute_sigma2(j.A);	\
	  sigma3_ ## A[im_r_ilinmom][ijack]=compute_sigma3(j.A)
	  
	  //compute sigma
	  sigma_time.start();
	  auto &j=jprops_inv[im_r];
	  COMPUTE_SIGMA(LO);
	  if(pars::use_QED)
	    {
	      COMPUTE_SIGMA(CR_CT);
	      COMPUTE_SIGMA(TM_CT);
	      COMPUTE_SIGMA(PH);
	    }
	  sigma_time.stop();
	  
#undef COMPUTE_SIGMA
	}
    }
  
  return *this;
}

void perens_t::average_r_sigma(perens_t &out) const
{
  cout<<"Averaging r for sigma"<<endl;
  
  if(nr==1)
    {
      cout<<"Skipping r average of Zq"<<endl;
      out=*this;
    }
  else
    for(auto &t : out.get_sigma_tasks({this}))
      {
	cout<<" "<<t.tag<<endl;
	
	const djvec_t &sigma=*t.in.front();
	djvec_t &sigma_rave=*t.out;
	
	for(size_t out_i=0;out_i<out.im_r_ilinmom_ind.max();out_i++)
	  {
	    const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(out_i);
	    vector<size_t> im_r_ilinmom_comp=out_im_r_ilinmom_comp;
	    
	    sigma_rave[out_i]=0.0;
	    for(size_t r=0;r<nr;r++)
	      {
		im_r_ilinmom_comp[1]=r;
		const size_t i=im_r_ilinmom_ind(im_r_ilinmom_comp);
		sigma_rave[out_i]+=sigma[i];
	      }
	    sigma_rave[out_i]/=nr;
	  }
      }
}

void perens_t::average_equiv_momenta_sigma(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const
{
  for(size_t i=0;i<out.im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(i);
      const size_t out_ilinmom_combo=out_im_r_ilinmom_comp[2];
      
      for(const auto &t : out.get_sigma_tasks({this}))
	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq : equiv_linmom_combos[out_ilinmom_combo])
	    {
	      vector<size_t> in_im_r_ilinmom_comp=out_im_r_ilinmom_comp;
	      in_im_r_ilinmom_comp[2]=ieq;
	      const size_t i=im_r_ilinmom_ind(in_im_r_ilinmom_comp);
	      
	      ave+=(*t.in.front())[i];
	    }
  	  ave/=equiv_linmom_combos[out_ilinmom_combo].size();
  	}
    }
}

void perens_t::val_chir_extrap_sigma(perens_t &out) const
{
  //slice m
  vector<double> x(nm);
  djvec_t y(nm);
  for(size_t im=0;im<nm;im++)
    if(pars::chir_extr_method==chir_extr::MQUARK) x[im]=am[im];
    else                                          x[im]=sqr(meson_mass[im_im_ind({im,im})].ave());
  
  for(auto &t : out.get_sigma_tasks({this}))
    {
      const djvec_t &sigma=*t.in.front();
      djvec_t &sigma_chir=*t.out;
      const string &tag=t.tag;
      
      for(size_t ilinmom=0;ilinmom<linmoms.size();ilinmom++)
	{
	  //open the plot file if needed
	  const string plot_path=dir_path+"/plots/chir_extr_"+tag+"_mom_"+to_string(ilinmom)+".xmg";
	  grace_file_t *plot=nullptr;
	  if(ilinmom%pars::print_each_mom==0) plot=new grace_file_t(plot_path);
	  
	  for(size_t r=0;r<nr;r++)
	    {
	      //slice m
	      djvec_t y(nm);
	      for(size_t im=0;im<nm;im++)
		y[im]=sigma[im_r_ilinmom_ind({im,r,ilinmom})];
	      
	      //fit, store and write the result
	      djvec_t coeffs=poly_fit(x,y,1);
	      sigma_chir[out.im_r_ilinmom_ind({0,r,ilinmom})]=coeffs[0];
	      if(plot!=nullptr)
		{
		  auto xminmax=minmax_element(x.begin(),x.end());
		  double xmax=*xminmax.second*1.1;
		  write_fit_plot(*plot,0,xmax,bind(poly_eval<djvec_t>,coeffs,_1),x,y);
		  plot->write_ave_err(0,coeffs[0].ave_err());
		}
	    }
	  
	  if(plot!=nullptr) delete plot;
	}
    }
}
