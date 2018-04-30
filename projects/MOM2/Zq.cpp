#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <MOM2/Zq.hpp>

#include <MOM2/perens.hpp>

void perens_t::compute_Zq()
{
#pragma omp parallel for
  for(size_t im_r_ilinmom=0;im_r_ilinmom<im_r_ilinmom_ind.max();im_r_ilinmom++)
    {
      const vector<size_t> comps=im_r_ilinmom_ind(im_r_ilinmom);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t im_r=im_r_ind({im,r});
      Zq[im_r_ilinmom]=sigma1_LO[im_r_ilinmom];
      
      if(pars::use_QED)
	{
	  Zq_QED[im_r_ilinmom]=
	    sigma1_PH[im_r_ilinmom]+
	    sigma1_CR_CT[im_r_ilinmom]*deltam_cr[im_r]+
	    sigma1_TM_CT[im_r_ilinmom]*deltam_tm[im_r];
	}
    }
}

vector<perens_t::task_t> perens_t::get_Zq_tasks(const vector<const perens_t*>& ens)
{
  vector<const djvec_t*> in_Zq,in_Zq_QED;
  for(auto &e : ens)
    {
      in_Zq.push_back(&e->Zq);
      if(pars::use_QED)
	in_Zq_QED.push_back(&e->Zq_QED);
    }
  
  vector<task_t> Zq_tasks={{&Zq,in_Zq,im_r_ilinmom_ind,"Zq",QCD_task}};
  if(pars::use_QED)
    Zq_tasks.push_back({&Zq_QED,in_Zq_QED,im_r_ilinmom_ind,"Zq_QED",QED_task});
  
  return Zq_tasks;
}

void perens_t::plot_Zq(const string &suffix)
{
  cout<<"Plotting all Zq of "<<dir_path<<" for suffix: \""<<suffix<<"\""<<endl;
  
  auto tasks=this->get_Zq_tasks();
  
  djvec_t Zq_QED_rel;
  if(pars::use_QED)
    {
      Zq_QED_rel=Zq_QED/Zq;
      
      tasks.push_back({&Zq_QED_rel,{},im_r_ilinmom_ind,"Zq_QED_rel",QED_task});
    }
  
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

void perens_t::average_equiv_momenta_Zq(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const
{
  for(size_t i=0;i<out.im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(i);
      const size_t out_ilinmom_combo=out_im_r_ilinmom_comp[2];
      
      for(const auto &t : out.get_Zq_tasks({this}))
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

void perens_t::val_chir_extrap_Zq(perens_t &out) const
{
  //slice m
  vector<double> x(nm);
  djvec_t y(nm);
  for(size_t im=0;im<nm;im++)
    if(pars::chir_extr_method==chir_extr::MQUARK) x[im]=am[im];
    else                                          x[im]=sqr(meson_mass[im_im_ind({im,im})].ave());
  
  for(auto &t : out.get_Zq_tasks({this}))
    {
      const djvec_t &Zq=*t.in.front();
      djvec_t &Zq_chir=*t.out;
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
		y[im]=Zq[im_r_ilinmom_ind({im,r,ilinmom})];
	      
	      //fit, store and write the result
	      djvec_t coeffs=poly_fit(x,y,1);
	      Zq_chir[out.im_r_ilinmom_ind({0,r,ilinmom})]=coeffs[0];
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
