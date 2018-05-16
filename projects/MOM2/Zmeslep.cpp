#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZMESLEP
 #include <MOM2/Zmeslep.hpp>

#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

using namespace meslep;

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
  cout<<"Computing Zmeslep"<<endl;
  
  for(size_t im_in=0;im_in<nm;im_in++)
    for(size_t r_in=0;r_in<nr;r_in++)
      for(size_t im_ou=0;im_ou<nm;im_ou++)
	for(size_t r_ou=0;r_ou<nr;r_ou++)
	  {
	    const size_t im_r_in=im_r_ind({im_in,r_in});
	    const size_t im_r_ou=im_r_ind({im_ou,r_ou});
	    
	    cout<<"deltam_cr[im_r_in="<<im_r_in<<", "<<im_r_ind.descr(im_r_in)<<"]: "<<deltam_cr[im_r_in].ave_err()<<endl;
	    cout<<"deltam_cr[im_r_ou="<<im_r_ou<<", "<<im_r_ind.descr(im_r_ou)<<"]: "<<deltam_cr[im_r_ou].ave_err()<<endl;
	    cout<<"deltam_tm[im_r_in="<<im_r_in<<", "<<im_r_ind.descr(im_r_in)<<"]: "<<deltam_tm[im_r_in].ave_err()<<endl;
	    cout<<"deltam_tm[im_r_ou="<<im_r_ou<<", "<<im_r_ind.descr(im_r_ou)<<"]: "<<deltam_tm[im_r_ou].ave_err()<<endl;
	    
#pragma omp parallel for
	    for(size_t imeslepmom=0;imeslepmom<meslepmoms().size();imeslepmom++)
	      {
		using namespace meslep;
		
		const size_t ilinmom_in=meslepmoms()[imeslepmom][1];
		const size_t ilinmom_ou=meslepmoms()[imeslepmom][2];
		auto s1_in=sigma_ins_getter(im_in,r_in,ilinmom_in,sigma::SIGMA1);
		auto s1_ou=sigma_ins_getter(im_ou,r_ou,ilinmom_ou,sigma::SIGMA1);
		
		djack_t sigma1_QED_in=
		  (s1_in(sigma::PH)+
		   s1_in(sigma::CR)*deltam_cr[im_r_in]+
		   s1_in(sigma::TM)*deltam_tm[im_r_in])*q_in*q_in;
		djack_t sigma1_QED_ou=
		  (s1_ou(sigma::PH)+
		   s1_ou(sigma::CR)*deltam_cr[im_r_ou]+
		   s1_ou(sigma::TM)*deltam_tm[im_r_ou])*q_ou*q_ou;
		
		for(size_t ijack=0;ijack<=njacks;ijack++)
		  {
		    using Zmeslep_t=Matrix<double,nZop,nZop>;
		    Zmeslep_t Gamma_meslep_combo;
		    Zmeslep_t Gamma_QED_meslep_combo;
		    
		    for(size_t iop=0;iop<nZop;iop++)
		      for(size_t iproj=0;iproj<nZop;iproj++)
			{
			  using namespace pr_meslep;
			  
			  auto ml=pr_meslep_ins_getter(im_ou,r_ou,im_in,r_in,iop,iproj,imeslepmom);
			  
			  //Gamma LO and correction
			  Gamma_meslep_combo(iop,iproj)=ml(LO)[ijack];
			  
			  //QED correction
			  Gamma_QED_meslep_combo(iop,iproj)=
			    ml(EX)[ijack]*q_in*q_ou+
			    
			    (ml(NA_IN)[ijack]*q_in+
			     ml(NA_OU)[ijack]*q_ou)*ql+
			    
			    (ml(PH_IN)[ijack]+
			     ml(CR_IN)[ijack]*deltam_cr[im_r_in][ijack]+
			     ml(TM_IN)[ijack]*deltam_tm[im_r_in][ijack])*q_in*q_in+
			    
			    (ml(PH_OU)[ijack] +
			     ml(CR_OU)[ijack]*deltam_cr[im_r_ou][ijack]+
			     ml(TM_OU)[ijack]*deltam_tm[im_r_ou][ijack])*q_ou*q_ou
			    ;
			}
		    
		    const Zmeslep_t Gamma_meslep_combo_inv=Gamma_meslep_combo.inverse();
		    
		    const double Zq_contr=sqrt(s1_ou(sigma::LO)[ijack]*s1_in(sigma::LO)[ijack]);
		    const double Zq_QED_rel_contr=(sigma1_QED_in[ijack]/s1_in(sigma::LO)[ijack]+
						   sigma1_QED_ou[ijack]/s1_ou(sigma::LO)[ijack])/2.0;
		    
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
}

void perens_t::plot_Zmeslep(const string &suffix)
{
  cout<<"Plotting all Zmeslep of "<<dir_path<<" for suffix: \""<<suffix<<"\""<<endl;
  
  for(const auto &t : this->get_Zmeslep_tasks())
    {
      //decript tuple
      const djvec_t &Z=*t.out;
      const string &tag=t.tag;
      
      for(size_t iop=0;iop<nbil;iop++)
	for(size_t iproj=0;iproj<nbil;iproj++)
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
