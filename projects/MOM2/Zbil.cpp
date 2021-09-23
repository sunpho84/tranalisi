#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_ZBIL
 #include <MOM2/Zbil.hpp>

#include <MOM2/analysis.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/timings.hpp>

vector<perens_t::task_t> perens_t::get_Zbil_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_Zbil,in_Zbil_QED_rel;
  for(auto &e : ens)
    {
      in_Zbil.push_back(&e->Zbil);
      if(pars::use_QED) in_Zbil_QED_rel.push_back(&e->Zbil_QED_rel);
    }
  vector<task_t> Zbil_tasks={{&Zbil,in_Zbil,im_r_im_r_ibil_ibilmom_ind,"Zbil",QCD_task}};
  if(pars::use_QED) Zbil_tasks.push_back({&Zbil_QED_rel,in_Zbil_QED_rel,im_r_im_r_ibil_ibilmom_ind,"Zbil"+QED_tag_suffix(),QED_task});
  
  return Zbil_tasks;
}

void perens_t::evolve_QED_Zbil_mixed_to_1_ov_a(perens_t& out) const
{
  cout<<"Evolving Zbil"<<endl;
  for(size_t im_r_im_r_ibil_ibilmom=0;im_r_im_r_ibil_ibilmom<im_r_im_r_ibil_ibilmom_ind.max();im_r_im_r_ibil_ibilmom++)
    {
      const vector<size_t> im_r_im_r_ibil_ibilmom_comps=im_r_im_r_ibil_ibilmom_ind(im_r_im_r_ibil_ibilmom);
      const size_t ibil=im_r_im_r_ibil_ibilmom_comps[4];
      const size_t ibil_mom=im_r_im_r_ibil_ibilmom_comps[5];
      const double a2p2=all_moms[bilmoms[ibil_mom][0]].p(L).norm2();
      const double p2=a2p2*sqr(ainv);
      
      //const double gamma[5]={-8.0,0.0,-8.0,0.0,-152.0/3};
      
      const double gamma_s0[nbil]  ={-6.0*CF,0.0,-6.0*CF,0.0,+2.0*CF};
      const double gamma_e0[nbil]  ={-6.0,   0.0,-6.0,   0.0,+2.0};
      const double gamma_se1[nbil] ={-8.0,   0.0,-8.0,   0.0,-152.0/3.0};
      
      out.Zbil_QED_rel[im_r_im_r_ibil_ibilmom]=
	Zbil_QED_rel[im_r_im_r_ibil_ibilmom]+
	evolve_QED_mixed_alpha(a2p2,gamma_s0[ibil],gamma_se1[ibil],gamma_e0[ibil],evol::evolution_Zbil_to_RIp(ibil_t_list[ibil],evol::NF4,/*ord*/3,ainv,p2));
    }
}

void perens_t::compute_Zbil(const bool also_QCD,const bool also_QED)
{
  cout<<"Computing Zbil"<<endl;
  
  if(also_QED and pars::use_QED) needs_to_read_assembled_QED_greenfunctions();
  
  for(size_t ibilmom=0;ibilmom<bilmoms.size();ibilmom++)
    for(size_t im_r_im_r_ibil=0;im_r_im_r_ibil<im_r_im_r_ibil_ind.max();im_r_im_r_ibil++)
      {
	const vector<size_t> im_r_im_r_ibil_comp=im_r_im_r_ibil_ind(im_r_im_r_ibil);
	const size_t im_ou=im_r_im_r_ibil_comp[0];
	const size_t r_ou=im_r_im_r_ibil_comp[1];
	const size_t im_in=im_r_im_r_ibil_comp[2];
	const size_t r_in=im_r_im_r_ibil_comp[3];
	const size_t ibil=im_r_im_r_ibil_comp[4];
	const size_t ilinmom_ou=bilmoms[ibilmom][1];
	const size_t ilinmom_in=bilmoms[ibilmom][2];
	
	const size_t im_r_im_r_ibil_ibilmom=im_r_im_r_ibil_ibilmom_ind(concat(im_r_im_r_ibil_comp,ibilmom));
	
	auto s1_ou=sigma_ins_getter(im_ou,r_ou,ilinmom_ou,sigma::SIGMA1);
	auto s1_in=sigma_ins_getter(im_in,r_in,ilinmom_in,sigma::SIGMA1);
	auto pr=pr_bil_ins_getter(im_ou,r_ou,im_in,r_in,ibil,ibilmom);
	
	if(also_QCD)
	  Zbil[im_r_im_r_ibil_ibilmom]=
	    sqrt(s1_ou(sigma::LO)*s1_in(sigma::LO))/pr(pr_bil::LO);
	
	if(also_QED and pars::use_QED)
	  {
	    Zbil_QED_rel[im_r_im_r_ibil_ibilmom]=
	      -pr(pr_bil::QED)/pr(pr_bil::LO);
	    
	    if(pars::include_Zq_corrections_in_bilinears)
	      Zbil_QED_rel[im_r_im_r_ibil_ibilmom]+=
		(s1_ou(sigma::QED)/s1_ou(sigma::LO)+
		 s1_in(sigma::QED)/s1_in(sigma::LO))/2.0;
	  }
      }
}

void perens_t::interpolate_Zbil_to_p2ref(perens_t &out) const
{
  cout<<"Interpolating to reference p2 Zbil"<<endl;
  
  //get ranges
  double a2p2=pars::p2ref/sqr(ainv);
  pair<double,double> a2p2minmax=get_a2p2tilde_range_bracketting(bilmoms,a2p2);
  const double p2min=a2p2minmax.first*sqr(ainv);
  const double p2max=a2p2minmax.second*sqr(ainv);
  cout<<"p2 range:   "<<p2min<<" - "<<p2max<<endl;
  
  //fit points
  vector<double> x(bilmoms.size());
  djvec_t y(bilmoms.size());

  for(auto &t : out.get_Zbil_tasks({this}))
    {
      const djvec_t &in=*t.in.front();
      djvec_t &out=*t.out;
      const string &tag=t.tag;
      
      if(fabs(in[0].ave()-1)<1e-7 and in[0].err()<1e-7)
	{
	  cout<<"Skipping interpolation for "<<tag<<endl;
	  out=in[0];
	}
      else
	for(size_t im_r_im_r_ibil=0;im_r_im_r_ibil<im_r_im_r_ibil_ind.max();im_r_im_r_ibil++)
	  {
	    const vector<size_t> comps=im_r_im_r_ibil_ind(im_r_im_r_ibil);
	    
	    //take physical units for plot
	    for(size_t imom=0;imom<bilmoms.size();imom++)
	      {
		vector<size_t> comps_with_mom=comps;
		comps_with_mom.push_back(imom);
		
		x[imom]=all_moms[bilmoms[imom][0]].p(L).tilde().norm2()*sqr(ainv);
		y[imom]=in[im_r_im_r_ibil_ibilmom_ind(comps_with_mom)];
	      }
	    
	    //fit and interpolate
	    const djvec_t coeffs=poly_fit(x,y,2,p2min,p2max);
	    out[im_r_im_r_ibil]=poly_eval(coeffs,pars::p2ref);
	    
	    //produce plot
	    const string path=dir_path+"/plots/interpolate_to_p2ref_"+tag+"_"+bil_tag[comps[4]]+
	      "_m"+to_string(comps[0])+"_r"+to_string(comps[1])+
	      "_m"+to_string(comps[2])+"_r"+to_string(comps[3])+".xmg";
	    grace_file_t plot(path);
	    write_fit_plot(plot,p2min,p2max,bind(poly_eval<djvec_t,double>,coeffs,_1),x,y);
	    plot.write_ave_err(pars::p2ref,out[im_r_im_r_ibil].ave_err());
	  }
    }
}

void perens_t::plot_Zbil(const string &suffix)
{
  cout<<"Plotting all Zbil of "<<dir_path<<" for suffix: \""<<suffix<<"\""<<endl;
  
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
  
  //plot also Zs/Zp
  if(pars::plot_Zs_fr_Zp)
    {
      grace_file_t out(dir_path+"/plots/Zs_fr_Zp"+(suffix!=""?("_"+suffix):string(""))+".xmg");
      grace_file_t out_QED;
      if(pars::use_QED) out_QED.open(dir_path+"/plots/Zs_fr_Zp_QED_rel"+(suffix!=""?("_"+suffix):string(""))+".xmg");
      
      //write mass by mass, only half of the combos
      for(size_t im1=0;im1<nm;im1++)
	for(size_t im2=im1;im2<nm;im2++)
	  for(size_t r=0;r<nr;r++)
	    {
	      out.new_data_set();
	      if(pars::use_QED) out_QED.new_data_set();
	      
	      for(size_t imom=0;imom<bilmoms.size();imom++)
		{
		  const double p2tilde=all_moms[bilmoms[imom][0]].p(L).tilde().norm2();
		  
		  const djack_t &ZS=Zbil[im_r_im_r_ibil_ibilmom_ind({im1,r,im2,r,iS,imom})];
		  const djack_t &ZP=Zbil[im_r_im_r_ibil_ibilmom_ind({im1,r,im2,r,iP,imom})];
		  const djack_t ZS_fr_ZP=ZS/ZP;
		  out.write_ave_err(p2tilde,ZS_fr_ZP.ave_err());
		  
		  if(pars::use_QED)
		    {
		      const djack_t &ZS_QED_rel=Zbil_QED_rel[im_r_im_r_ibil_ibilmom_ind({im1,r,im2,r,iS,imom})];
		      const djack_t &ZP_QED_rel=Zbil_QED_rel[im_r_im_r_ibil_ibilmom_ind({im1,r,im2,r,iP,imom})];
		      const djack_t ZS_fr_ZP_QED_rel=ZS_QED_rel-ZP_QED_rel;
		      out_QED.write_ave_err(p2tilde,ZS_fr_ZP_QED_rel.ave_err());
		    }
		}
	    }
    }
}
