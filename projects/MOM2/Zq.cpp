#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <MOM2/Zq.hpp>

#include <MOM2/perens.hpp>

double perens_t::compute_Zq(const qprop_t &prop_inv,const size_t glb_mom)
{
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const qprop_t pslash=qua_slash(ptilde);
  
  const double Zq=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}

djack_t perens_t::compute_Zq(const jqprop_t &jprop_inv,const size_t glb_mom)
{
  djack_t Zq;
  
  const p_t ptilde=all_moms[glb_mom].p(L).tilde();
  const double pt2=ptilde.norm2();
  const qprop_t pslash=qua_slash(ptilde);
  
  for(size_t ijack=0;ijack<=njacks;ijack++)
    Zq[ijack]=(jprop_inv[ijack]*pslash).trace().imag()/(12.0*pt2*V);
  
  return Zq;
}

void perens_t::plot_Zq(const string &suffix)
{
  for(auto &t : get_Zq_tasks(*this))
    {
      const djvec_t &Z=*t.in;
      const string &tag=t.tag;
      
      grace_file_t out(dir_path+"/plots/"+tag+(suffix!=""?("_"+suffix):string(""))+".xmg");
      
      for(size_t im=0;im<nm;im++)
  	for(size_t r=0;r<nr;r++)
	  {
	    out.new_data_set();
	    for(size_t imom=0;imom<linmoms.size();imom++)
	      {
		const double p2hat=all_moms[linmoms[imom][0]].p(L).tilde().norm2();
		out.write_ave_err(p2hat,Z[im_r_ilinmom_ind({im,r,imom})].ave_err());
	      }
	  }
    }
}

void perens_t::average_r_Zq(perens_t &out) const
{
  for(auto &t : get_Zq_tasks(out))
    {
      const djvec_t &Zq=*t.in;
      djvec_t &Zq_rave=*t.out;
      
      for(size_t out_i=0;out_i<out.im_r_ilinmom_ind.max();out_i++)
	{
	  const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(out_i);
	  vector<size_t> im_r_ilinmom_comp=out_im_r_ilinmom_comp;
	  
	  Zq_rave[out_i]=0.0;
	  for(size_t r=0;r<nr;r++)
	    {
	      im_r_ilinmom_comp[1]=r;
	      const size_t i=im_r_ilinmom_ind(im_r_ilinmom_comp);
	      Zq_rave[out_i]+=Zq[i];
	    }
	  Zq_rave[out_i]/=nr;
	}
    }
}

void perens_t::average_equiv_momenta_Zq(perens_t &out,const vector<vector<size_t>> &equiv_linmom_combos) const
{
  for(size_t i=0;i<out.im_r_ilinmom_ind.max();i++)
    {
      const vector<size_t> out_im_r_ilinmom_comp=out.im_r_ilinmom_ind(i);
      const size_t out_ilinmom_combo=out_im_r_ilinmom_comp[2];
      
      for(const auto &t : get_Zq_tasks(out))
	{
  	  djack_t &ave=(*t.out)[i];
  	  ave=0.0;
  	  for(const size_t ieq : equiv_linmom_combos[out_ilinmom_combo])
	    {
	      vector<size_t> in_im_r_ilinmom_comp=out_im_r_ilinmom_comp;
	      in_im_r_ilinmom_comp[2]=ieq;
	      const size_t i=im_r_ilinmom_ind(in_im_r_ilinmom_comp);
	      
	      ave+=(*t.in)[i];
	    }
  	  ave/=equiv_linmom_combos[out_ilinmom_combo].size();
  	}
    }
}

void perens_t::val_chir_extrap_Zq(perens_t &out) const
{
  for(auto &t : get_Zq_tasks(out))
    {
      const djvec_t &Zq=*t.in;
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
	      for(size_t im=0;im<nm;im++) y[im]=Zq[im_r_ilinmom_ind({im,r,ilinmom})];
	      
	      //fit, store and write the result
	      djvec_t coeffs=poly_fit(am,y,1,am_min(),am_max());
	      Zq_chir[out.im_r_ilinmom_ind({0,r,ilinmom})]=coeffs[0];
	      if(plot!=nullptr)
		{
		  write_fit_plot(*plot,0,am_max(),bind(poly_eval<djvec_t>,coeffs,_1),am,y);
		  plot->write_ave_err(0,coeffs[0].ave_err());
		}
	    }
	  
	  if(plot!=nullptr) delete plot;
	}
    }
}
