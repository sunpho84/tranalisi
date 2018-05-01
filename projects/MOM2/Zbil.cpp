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

void perens_t::compute_Zbil()
{
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(size_t im_r_im_r_ibil=0;im_r_im_r_ibil<im_r_im_r_ibil_ind.max();im_r_im_r_ibil++)
      {
	CRASH("");
	
	// const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
	// const vector<size_t> im_r1_comp=subset(im_r_im_r_ibil_comp,0,2);
	// const vector<size_t> im_r2_comp=subset(im_r_im_r_ibil_comp,2,4);
	// const size_t ibil=im_r_im_r_ibil_comp[4];
	// const size_t ilinmom1=bilmoms[ibilmom][1];
	// const size_t ilinmom2=bilmoms[ibilmom][2];
	// const size_t im_r1_ilinmom1=im_r_ilinmom_ind(concat(im_r1_comp,ilinmom1));
	// const size_t im_r2_ilinmom2=im_r_ilinmom_ind(concat(im_r2_comp,ilinmom2));
	
	// const size_t im_r_im_r_ibil_ibilmom=im_r_im_r_ibil_ibilmom_ind(concat(im_r1_comp,im_r2_comp,vector<size_t>({ibil,ibilmom})));
	
	// bil[im_r_im_r_ibil_ibilmom]=
	//   sqrt(Zq_sig1[im_r1_ilinmom1]*Zq_sig1[im_r2_ilinmom2])/pr_bil[im_r_im_r_ibil_ibilmom];
	
	// if(pars::use_QED)
	//     bil_QED_rel[im_r_im_r_ibil_ibilmom]=
	//       -pr_bil_QED[im_r_im_r_ibil_ibilmom]/pr_bil[im_r_im_r_ibil_ibilmom]
	//       +(Zq_sig1_QED[im_r1_ilinmom1]/Zq_sig1[im_r1_ilinmom1]+
	// 	Zq_sig1_QED[im_r2_ilinmom2]/Zq_sig1[im_r2_ilinmom2])/2.0;
      }
}

void perens_t::plot_Zbil(const string &suffix)
{
  for(const auto &t : this->get_Zbil_tasks())
    {
      //decript tuple
      const djvec_t &Z=*t.out;
      const string &tag=t.tag;
      
      for(size_t ibil=0;ibil<nbil;ibil++)
  	{
	  grace_file_t out(dir_path+"/plots/"+tag+"_"+bil_tag[ibil]+(suffix!=""?("_"+suffix):string(""))+".xmg");
	  
	  //write mass by mass, only half of the combos
  	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=im1;im2<nm;im2++)
	      for(size_t r=0;r<nr;r++)
		{
		  out.new_data_set();
		  
		  for(size_t imom=0;imom<bilmoms.size();imom++)
		    {
		      const double p2tilde=all_moms[bilmoms[imom][0]].p(L).tilde().norm2();
		      out.write_ave_err(p2tilde,Z[im_r_im_r_ibil_ibilmom_ind({im1,r,im2,r,ibil,imom})].ave_err());
		    }
		}
	}
    }
}
