#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <contractions.hpp>
#include <deltam_cr.hpp>
#include <corrections.hpp>
#include <evolutions.hpp>
#include <geometry.hpp>
#include <ingredients.hpp>
#include <prop.hpp>
#include <read.hpp>
#include <sig3.hpp>
#include <timings.hpp>
#include <types.hpp>
#include <Zbil.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

//! write a given Z
void write_Z(const string &name,const djvec_t &Z,const vector<double> &pt2)
{
  grace_file_t outf("plots/"+name+".xmg");
  outf.write_vec_ave_err(pt2,Z.ave_err());
}

//! linearly fit a given Z
djvec_t linfit_Z(const djvec_t &Z,const string &name,double band_val=0)
{
  vector<double> pt2=get_indep_pt2();
  double p2max=*max_element(pt2.begin(),pt2.end());
  djvec_t pars=poly_fit(pt2,Z,1,1.0,2.0);
  
  grace_file_t outf("plots/"+name+".xmg");
  outf.write_vec_ave_err(pt2,Z.ave_err());
  outf.write_polygon(bind(poly_eval<djvec_t>,pars,_1),0,p2max);
  if(band_val!=0) outf.write_line([band_val](double x){return band_val;},0,p2max);
  
  return pars;
}

int main(int narg,char **arg)
{
  //read input file
  string input_path="input.txt";
  if(narg>=2) input_path=arg[1];
  read_input(input_path);

  get_deltam_cr();
  
  const string ingredients_path="ingredients.dat";
  ingredients_t ing;
  if(file_exists(ingredients_path)) ing.bin_read(ingredients_path);
  else
    {
      ing.ri_mom();
      ing.bin_write(ingredients_path);
    }
  
  //read or write?
  
  //Subtracted Zq, with and without EM, all moms, averaged r
  djvec_t Zq_chir_allmoms(imoms.size());
  djvec_t Zq_chir_allmoms_sub(imoms.size());
  djvec_t Zq_chir_allmoms_sub_evolved(imoms.size());
  djvec_t Zq_sig1_chir_allmoms(imoms.size());
  djvec_t Zq_sig1_chir_allmoms_sub(imoms.size());
  djvec_t Zq_sig1_chir_allmoms_sub_evolved(imoms.size());
  djvec_t Zq_sig1_EM_chir_allmoms(imoms.size());
  djvec_t Zq_sig1_EM_chir_allmoms_sub(imoms.size());
  djvec_t Zq_sig1_EM_chir_allmoms_sub_evolved(imoms.size());
  
  //! list of task to chirally extrapolate Zq
  vector<tuple<djvec_t*,djvec_t*,string>> Zq_chirextr_tasks{
    {&ing.Zq_allmoms,&Zq_chir_allmoms,string("Zq")},
    {&ing.Zq_sig1_allmoms,&Zq_sig1_chir_allmoms,"Zq_sig1"}};
  if(use_QED) Zq_chirextr_tasks.push_back(make_tuple(&ing.Zq_sig1_EM_allmoms,&Zq_sig1_EM_chir_allmoms,"Zq_sig1_EM"));
  
  djvec_t Zbil_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms(iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms_sub(iZbil_imom_ind.max());
  djvec_t Zbil_chir_allmoms_sub_evolved(iZbil_imom_ind.max());
  djvec_t Zbil_QED_allmoms(im_r_im_r_iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms(iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms_sub(iZbil_imom_ind.max());
  djvec_t Zbil_QED_chir_allmoms_sub_evolved(iZbil_imom_ind.max());
  
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      //extrapolate to chiral limit Zq
      for(auto & p : Zq_chirextr_tasks)
	{
	  djvec_t y(nm);
	  const djvec_t &Zq=(*get<0>(p));
	  djvec_t &Zq_chir=(*get<1>(p));
	  const string &tag=get<2>(p);
	  
	  //slice m
	  double am_max=*max_element(am.begin(),am.end())*1.1;
	  for(size_t im=0;im<nm;im++)
	    {
	      //averages r if both asked
	      y[im]=0.0;
	      for(size_t r=0;r<nr;r++)
		{
		  size_t i=im_r_imom_ind({im,r,imom});
		  y[im]+=Zq[i]/nr;
		  cout<<tag<<"["<<i<<"=(im"<<im<<"r"<<r<<"imom"<<imom<<")], mom "<<imoms[imom].p(L).norm2()<<": "<<Zq[i].ave_err()<<endl;
		}
	      cout<<tag<<"[im"<<im<<"imom"<<imom<<"]: "<<y[im][0]<<" "<<y[im].err()<<endl;
	    }
	  
	  //fit and write the result
	  djvec_t coeffs=poly_fit(am,y,1,am_min,am_max);
	  if(imom%print_each_mom==0)
	    {
	      grace_file_t plot("plots/chir_extr_"+tag+"_mom_"+to_string(imom)+".xmg");
	      write_fit_plot(plot,0,am_max,bind(poly_eval<djvec_t>,coeffs,_1),am,y);
	      plot.write_ave_err(0,coeffs[0].ave_err());
	    }
	  //extrapolated value
	  Zq_chir[imom]=coeffs[0];
	}
      
      //subtract Zq
      imom_t mom=imoms[imom];
      double sub=g2tilde*sig1_a2(act,gf::LANDAU,group::SU3,mom,L);
      Zq_chir_allmoms_sub[imom]=Zq_chir_allmoms[imom]-sub;
      Zq_sig1_chir_allmoms_sub[imom]=Zq_sig1_chir_allmoms[imom]-sub;
      if(use_QED)
	{
	  double sub_EM=-1.0*/*coupling gets -1 due to definition of expansion*/sig1_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L);
	  Zq_sig1_EM_chir_allmoms_sub[imom]=Zq_sig1_EM_chir_allmoms[imom]-sub_EM*
	    Zq_sig1_chir_allmoms[imom]; //factorization hypotesis
	}
      
      //evolver
      double p2=mom.p(L).norm2();
      double evolver_Zq=evolution_Zq_to_RIp(Nf,ord,ainv,p2);
      cout<<"EvolverZq["<<p2<<"]="<<evolver_Zq<<endl;
      
      Zq_chir_allmoms_sub_evolved[imom]=Zq_chir_allmoms_sub[imom]/evolver_Zq;
      Zq_sig1_chir_allmoms_sub_evolved[imom]=Zq_sig1_chir_allmoms_sub[imom]/evolver_Zq;
      Zq_sig1_EM_chir_allmoms_sub_evolved[imom]=Zq_sig1_EM_chir_allmoms_sub[imom]/evolver_Zq;
      
      //chiral extrapolate
      djvec_t pr_bil_chir_mom(nZbil);
      djvec_t pr_bil_QED_chir_mom(nZbil);
      vector<tuple<djvec_t*,djvec_t*,string>>
	pr_bil_chirextr_tasks{{&ing.pr_bil_mom,&pr_bil_chir_mom,string("pr_bil")}};
      if(use_QED) pr_bil_chirextr_tasks.push_back(make_tuple(&ing.pr_bil_QED_mom,&pr_bil_QED_chir_mom,string("pr_bil_QED")));
      
      for(auto & p : pr_bil_chirextr_tasks)
	for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	  {
	    const djvec_t &pr=(*get<0>(p));
	    djvec_t &pr_chir=(*get<1>(p));
	    const string &tag=get<2>(p);
	    
	    //check if we need to subtract the pole
	    const bool sub_pole=(iZbil==iZS or iZbil==iZP);
	    
	    //slice m and fit it
	    djvec_t y(nm*(nm+1)/2),y_plot(nm*(nm+1)/2);
	    vector<double> x(nm*(nm+1)/2);
	    int i=0;
	    for(size_t im1=0;im1<nm;im1++)
	      for(size_t im2=im1;im2<nm;im2++)
		{
		  //compute mass sum
		  x[i]=am[im1]+am[im2];
		  //compute y and y_plot
		  y_plot[i]=0.0;
		  for(size_t r=0;r<nr;r++) y_plot[i]+=pr[im_r_im_r_iZbil_imom_ind({im1,r,im2,r,iZbil,imom})]/nr;
		  
		  if(sub_pole) y[i]=x[i]*y_plot[i];
		  else         y[i]=y_plot[i];
		  //increment
		  i++;
		}
	    
	    //fit and store extrapolated value
	    djvec_t coeffs=poly_fit(x,y,(sub_pole?2:1),2.0*am_min,2.0*am_max);
	    pr_chir[iZbil]=coeffs[sub_pole?1:0];
	    
	    //plot
	    if(imom%print_each_mom==0)
	      {
		grace_file_t plot("plots/chir_extr_"+tag+"_"+Zbil_tag[iZbil]+"_mom_"+to_string(imom)+".xmg");
		write_fit_plot(plot,2*am_min,2*am_max,[&coeffs,sub_pole](double x)->djack_t{return poly_eval<djvec_t>(coeffs,x)/(sub_pole?x:1);},x,y_plot);
		plot.write_ave_err(0.0,pr_chir[iZbil].ave_err());
	      }
	  }
      
      //compute subtractions
      djvec_t pr_bil_mom_correction(nZbil);
      djvec_t pr_bil_QED_mom_correction(nZbil);
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  pr_bil_mom_correction[iZbil]=g2tilde*pr_bil_a2(act,gf::LANDAU,group::SU3,mom,L,iZbil);
	  if(use_QED) pr_bil_QED_mom_correction[iZbil]=1.0*pr_bil_a2(gaz::PLAQ,gf::FEYNMAN,group::U1,mom,L,iZbil)*pr_bil_chir_mom[iZbil]; //factorization hypotesis
	}
      
      //subtract cut-off effects
      djvec_t pr_bil_chir_mom_sub(nZbil);
      djvec_t pr_bil_QED_chir_mom_sub(nZbil);
      vector<tuple<djvec_t*,djvec_t*,djvec_t*,string>>
	pr_bil_sub_tasks{{&pr_bil_chir_mom,&pr_bil_chir_mom_sub,&pr_bil_mom_correction,string("pr_bil")}};
      if(use_QED) pr_bil_sub_tasks.push_back({&pr_bil_QED_chir_mom,&pr_bil_QED_chir_mom_sub,&pr_bil_QED_mom_correction,string("pr_bil_QED")});
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	for(auto & p : pr_bil_sub_tasks)
	  {
	    djvec_t &pr_chir=(*get<0>(p));
	    djvec_t &pr_chir_sub=(*get<1>(p));
	    djvec_t &pr_corr=(*get<2>(p));
	    const string &tag=get<3>(p);
	    
	    //subtract from bilinear
	    pr_chir_sub[iZbil]=pr_chir[iZbil]-pr_corr[iZbil];
	    
	    //plot
	    if(imom%print_each_mom==0)
	      {
		grace_file_t plot("plots/sub_"+tag+"_"+Zbil_tag[iZbil]+"_mom_"+to_string(imom)+".xmg");
		plot.write_ave_err(0.0,pr_chir[iZbil].ave_err());
		plot.new_data_set();
		plot.write_ave_err(0.0,pr_chir_sub[iZbil].ave_err());
	      }
	}
      
      //build Z
      for(size_t im_r_im_r_iZbil=0;im_r_im_r_iZbil<im_r_im_r_iZbil_ind.max();im_r_im_r_iZbil++)
	{
	  const vector<size_t> im_r_im_r_iZbil_comp=im_r_im_r_iZbil_ind(im_r_im_r_iZbil);
	  const vector<size_t> im_r1_comp=subset(im_r_im_r_iZbil_comp,0,2);
	  const vector<size_t> im_r2_comp=subset(im_r_im_r_iZbil_comp,2,4);
	  const size_t iZbil=im_r_im_r_iZbil_comp[4];
	  const size_t im_r1_imom=im_r_imom_ind(concat(im_r1_comp,imom));
	  const size_t im_r2_imom=im_r_imom_ind(concat(im_r2_comp,imom));
	  
	  const size_t im_r_im_r_iZbil_imom=im_r_im_r_iZbil_imom_ind(concat(im_r1_comp,im_r2_comp,vector<size_t>({iZbil,imom})));
	  
	  Zbil_allmoms[im_r_im_r_iZbil_imom]=
	    sqrt(ing.Zq_sig1_allmoms[im_r1_imom]*ing.Zq_sig1_allmoms[im_r2_imom])/ing.pr_bil_mom[im_r_im_r_iZbil_imom];
	  
	  if(use_QED)
	    {
	      Zbil_QED_allmoms[im_r_im_r_iZbil_imom]=
		ing.pr_bil_QED_mom[im_r_im_r_iZbil]/ing.pr_bil_mom[im_r_im_r_iZbil]+
		(ing.Zq_sig1_EM_allmoms[im_r1_imom]/ing.Zq_sig1_allmoms[im_r1_imom]+ing.Zq_sig1_EM_allmoms[im_r2_imom]/
		 ing.Zq_sig1_allmoms[im_r2_imom])/2.0;
	    }
	}
      
      //builds Z in the chiral limit and evolves it
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  const size_t iZbil_imom=iZbil_imom_ind({iZbil,imom});
	  
	  double p2=mom.p(L).norm2();
	  const double evolver_Zbil=evolution_Zbil_to_RIp(iZbil_t_list[iZbil],Nf,ord,ainv,p2);
	  cout<<"EvolverZ"<<Zbil_tag[iZbil]<<"["<<p2<<"]="<<evolver_Zbil<<endl;
	  Zbil_chir_allmoms[iZbil_imom]=sqrt(Zq_sig1_chir_allmoms[imom]*Zq_sig1_chir_allmoms[imom])/pr_bil_chir_mom[iZbil];
	  Zbil_chir_allmoms_sub[iZbil_imom]=sqrt(Zq_sig1_chir_allmoms_sub[imom]*Zq_sig1_chir_allmoms_sub[imom])/pr_bil_chir_mom_sub[iZbil];
	  Zbil_chir_allmoms_sub_evolved[iZbil_imom]=Zbil_chir_allmoms_sub[iZbil_imom]/evolver_Zbil;
	  
	  if(use_QED)
	    {
	      Zbil_QED_chir_allmoms[iZbil_imom]=
		pr_bil_QED_chir_mom[iZbil]/pr_bil_chir_mom[iZbil]+
		    (Zq_sig1_EM_chir_allmoms[imom]/Zq_sig1_chir_allmoms[imom]+Zq_sig1_EM_chir_allmoms[imom]/Zq_sig1_chir_allmoms[imom])/2.0;
	      Zbil_QED_chir_allmoms_sub[iZbil_imom]=
		pr_bil_QED_chir_mom_sub[iZbil]/pr_bil_chir_mom_sub[iZbil]+
		    (Zq_sig1_EM_chir_allmoms_sub[imom]/Zq_sig1_chir_allmoms_sub[imom]+Zq_sig1_EM_chir_allmoms_sub[imom]/Zq_sig1_chir_allmoms_sub[imom])/2.0;
	      Zbil_QED_chir_allmoms_sub_evolved[iZbil_imom]=
		Zbil_QED_chir_allmoms_sub[iZbil_imom]/evolver_Zbil;
	    }
	}
    }
  
  const index_t iZbil_indep_imom_ind({{"Zbil",nZbil},{"indep_mom",equiv_imoms.size()}});
  const index_t im_r_im_r_iZbil_indep_imom_ind=im_r_im_r_iZbil_ind*index_t({{"indep_mom",equiv_imoms.size()}});
  
  //average equiv moms
  djvec_t Zq=average_equiv_moms(ing.Zq_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zq_sig1=average_equiv_moms(ing.Zq_sig1_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zq_sig1_EM=average_equiv_moms(ing.Zq_sig1_EM_allmoms,im_r_indep_imom_ind,im_r_imom_ind);
  djvec_t Zbil=average_equiv_moms(Zbil_allmoms,im_r_im_r_iZbil_indep_imom_ind,im_r_im_r_iZbil_imom_ind);
  djvec_t Zbil_QED=use_QED?average_equiv_moms(Zbil_QED_allmoms,im_r_im_r_iZbil_indep_imom_ind,im_r_im_r_iZbil_imom_ind):djvec_t();
  
  //chirally extrapolated ones
  djvec_t Zq_chir=average_equiv_moms(Zq_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_chir_sub=average_equiv_moms(Zq_chir_allmoms_sub,indep_imom_ind,imom_ind);
  djvec_t Zq_chir_sub_evolved=average_equiv_moms(Zq_chir_allmoms_sub_evolved,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_chir=average_equiv_moms(Zq_sig1_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_EM_chir=average_equiv_moms(Zq_sig1_EM_chir_allmoms,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_chir_sub=average_equiv_moms(Zq_sig1_chir_allmoms_sub,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_chir_sub_evolved=average_equiv_moms(Zq_sig1_chir_allmoms_sub_evolved,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_EM_chir_sub=average_equiv_moms(Zq_sig1_EM_chir_allmoms_sub,indep_imom_ind,imom_ind);
  djvec_t Zq_sig1_EM_chir_sub_evolved=average_equiv_moms(Zq_sig1_EM_chir_allmoms_sub_evolved,indep_imom_ind,imom_ind);
  
  djvec_t Zbil_chir=average_equiv_moms(Zbil_chir_allmoms,iZbil_indep_imom_ind,iZbil_imom_ind);
  djvec_t Zbil_chir_sub=average_equiv_moms(Zbil_chir_allmoms_sub,iZbil_indep_imom_ind,iZbil_imom_ind);
  djvec_t Zbil_chir_sub_evolved=average_equiv_moms(Zbil_chir_allmoms_sub_evolved,iZbil_indep_imom_ind,iZbil_imom_ind);
  djvec_t Zbil_QED_chir=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms,iZbil_indep_imom_ind,iZbil_imom_ind):djvec_t();
  djvec_t Zbil_QED_chir_sub=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms_sub,iZbil_indep_imom_ind,iZbil_imom_ind):djvec_t();
  djvec_t Zbil_QED_chir_sub_evolved=use_QED?average_equiv_moms(Zbil_QED_chir_allmoms_sub_evolved,iZbil_indep_imom_ind,iZbil_imom_ind):djvec_t();
  
  using Z_plot_task_t=tuple<djvec_t*,djvec_t*,djvec_t*,djvec_t*,string>;
  
  //! list of task to plot chiral extrapolation Zq
  vector<Z_plot_task_t> Zq_plot_tasks{
    {&Zq,&Zq_chir,&Zq_chir_sub,&Zq_chir_sub_evolved,string("Zq")},
    {&Zq_sig1,&Zq_sig1_chir,&Zq_sig1_chir_sub,&Zq_sig1_chir_sub_evolved,"Zq_sig1"}};
  if(use_QED) Zq_plot_tasks.push_back(make_tuple(&Zq_sig1_EM,&Zq_sig1_EM_chir,&Zq_sig1_EM_chir_sub,&Zq_sig1_EM_chir_sub_evolved,"Zq_sig1_EM"));
  
  for(auto &p : Zq_plot_tasks)
    {
      //decript the tuple
      const djvec_t &Zq=(*get<0>(p));
      const djvec_t &Zq_chir=(*get<1>(p));
      const djvec_t &Zq_chir_sub=(*get<2>(p));
      const djvec_t &Zq_chir_sub_evolved=(*get<3>(p));
      const string &tag=get<4>(p);
      
      //loop over all r
      for(size_t r=0;r<nr;r++)
	{
	  //open the file
	  grace_file_t out("plots/"+tag+"_r_"+to_string(r)+".xmg");
	  out.new_data_set();
	  
	  //m
	  for(size_t im=0;im<nm+3;im++)
	    {
	      for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		{
		  djack_t y;
		  if(im<nm) y=Zq[im_r_indep_imom_ind({im,r,indep_imom})];
		  if(im==nm+0) y=Zq_chir[indep_imom];
		  if(im==nm+1) y=Zq_chir_sub[indep_imom];
		  if(im==nm+2) y=Zq_chir_sub_evolved[indep_imom];
		  
		  out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),y.ave_err());
		}
	      out.new_data_set();
	    }
	}
    }
  
  //! list of task to print the chiral extrapolate bilinears
  vector<Z_plot_task_t> Zbil_tasks{{&Zbil,&Zbil_chir,&Zbil_chir_sub,&Zbil_chir_sub_evolved,string("Zbil")}};
  if(use_QED) Zbil_tasks.push_back(make_tuple(&Zbil_QED,&Zbil_QED_chir,&Zbil_QED_chir_sub,&Zbil_QED_chir_sub_evolved,"Zbil_EM"));
  
  for(auto &p : Zbil_tasks)
    {
      //decript tuple
      const djvec_t &Z=(*get<0>(p));
      const djvec_t &Z_chir=(*get<1>(p));
      const djvec_t &Z_chir_sub=(*get<2>(p));
      const djvec_t &Z_chir_sub_evolved=(*get<3>(p));
      const string &tag=get<4>(p);
      
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	{
	  //open the file
	  grace_file_t out("plots/"+tag+"_Z"+Zbil_tag[iZbil]+".xmg");
	  out.new_data_set();
	  
	  //write mass by mass, only half of the combos
	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=im1;im2<nm;im2++)
	      {
		for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		  {
		    djack_t y;
		    y=0.0;
		    for(size_t r=0;r<nr;r++) y+=Z[im_r_im_r_iZbil_indep_imom_ind({im1,r,im2,r,iZbil,indep_imom})]/nr;
		    
		    out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),y.ave_err());
		  }
		out.new_data_set();
	      }
	  
	  //write chiral extrap and subtracted
	  for(auto &Ztag : vector<tuple<const djvec_t*,string>>{{&Z_chir,"chir"},{&Z_chir_sub,"sub"},{&Z_chir_sub_evolved,"evo"}})
	    {
	      out.set_legend(get<1>(Ztag));
	      for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
		out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),(*get<0>(Ztag))[iZbil_indep_imom_ind({iZbil,indep_imom})].ave_err());
	      out.new_data_set();
	    }
	}
    }
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
