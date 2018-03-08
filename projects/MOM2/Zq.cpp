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

