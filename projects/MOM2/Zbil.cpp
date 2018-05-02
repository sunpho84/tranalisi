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
	const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
	const vector<size_t> im_r_ou_comp=subset(im_r_im_r_ibil_comp,0,2);
	const vector<size_t> im_r_in_comp=subset(im_r_im_r_ibil_comp,2,4);
	const size_t im_r_ou=im_r_ind(im_r_ou_comp);
	const size_t im_r_in=im_r_ind(im_r_in_comp);
	const size_t ilinmom_ou=bilmoms[ibilmom][1];
	const size_t ilinmom_in=bilmoms[ibilmom][2];
	const size_t im_r_ilinmom_ou=im_r_ilinmom_ind(concat(im_r_ou_comp,ilinmom_ou));
	const size_t im_r_ilinmom_in=im_r_ilinmom_ind(concat(im_r_in_comp,ilinmom_in));
	
	const size_t im_r_im_r_ibil_ibilmom=im_r_im_r_ibil_ibilmom_ind(concat(im_r_im_r_ibil_comp,ibilmom));
	
	Zbil[im_r_im_r_ibil_ibilmom]=
	  sqrt(sigma1_LO[im_r_ilinmom_ou]*sigma1_LO[im_r_ilinmom_in])/pr_bil_LO[im_r_im_r_ibil_ibilmom];
	
	if(pars::use_QED)
	  {
	    if(not deltam_computed) CRASH("Needs to have computed deltam");
	    
	    djack_t pr_bil_QED=
	      pr_bil_PH[im_r_im_r_ibil_ibilmom]+
	      pr_bil_CR_CT_OU[im_r_im_r_ibil_ibilmom]*deltam_cr[im_r_ou]+
	      pr_bil_CR_CT_IN[im_r_im_r_ibil_ibilmom]*deltam_cr[im_r_in]+
	      pr_bil_TM_CT_OU[im_r_im_r_ibil_ibilmom]*deltam_tm[im_r_ou]+
	      pr_bil_TM_CT_IN[im_r_im_r_ibil_ibilmom]*deltam_tm[im_r_in];
	    djack_t sigma1_QED_ou=
	      sigma1_PH[im_r_ilinmom_ou]+
	      sigma1_CR_CT[im_r_ilinmom_ou]*deltam_cr[im_r_ou]+
	      sigma1_TM_CT[im_r_ilinmom_ou]*deltam_tm[im_r_ou];
	    djack_t sigma1_QED_in=
	      sigma1_PH[im_r_ilinmom_in]+
	      sigma1_CR_CT[im_r_ilinmom_in]*deltam_cr[im_r_in]+
	      sigma1_TM_CT[im_r_ilinmom_in]*deltam_tm[im_r_in];
	    
	    Zbil_QED_rel[im_r_im_r_ibil_ibilmom]=
	      -pr_bil_QED/pr_bil_LO[im_r_im_r_ibil_ibilmom]*0.0
	      +(sigma1_QED_ou/sigma1_LO[im_r_ilinmom_ou]+
		sigma1_QED_in/sigma1_LO[im_r_ilinmom_in])/2.0;
	  }
      }
}

vector<perens_t::task_t> perens_t::get_Zbil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_Zbil,in_Zbil_QED_rel;
  for(auto &e : ens)
    {
      in_Zbil.push_back(&e->Zbil);
      if(pars::use_QED) in_Zbil_QED_rel.push_back(&e->Zbil_QED_rel);
    }
  vector<task_t> Zbil_tasks={{&Zbil,in_Zbil,im_r_im_r_ibil_ibilmom_ind,"Zbil",QCD_task}};
  if(pars::use_QED) Zbil_tasks.push_back({&Zbil_QED_rel,in_Zbil_QED_rel,im_r_im_r_ibil_ibilmom_ind,"Zbil_QED_rel",QED_task});
  
  return Zbil_tasks;
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
