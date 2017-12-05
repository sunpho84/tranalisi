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
// djvec_t linfit_Z(const djvec_t &Z,const string &name,double band_val=0)
// {
//   vector<double> pt2=get_indep_pt2();
//   double p2max=*max_element(pt2.begin(),pt2.end());
//   djvec_t pars=poly_fit(pt2,Z,1,1.0,2.0);
  
//   grace_file_t outf("plots/"+name+".xmg");
//   outf.write_vec_ave_err(pt2,Z.ave_err());
//   outf.write_polygon(bind(poly_eval<djvec_t>,pars,_1),0,p2max);
//   if(band_val!=0) outf.write_line([band_val](double x){return band_val;},0,p2max);
  
//   return pars;
// }

int main(int narg,char **arg)
{
  //read input file
  string input_path="input.txt";
  if(narg>=2) input_path=arg[1];
  read_input(input_path);
  
  get_deltam_cr();
  
  ingredients_t ing;
  ing.create_from_scratch();
  
  ing.plot_Z();
  
  ingredients_t chir=ing.chir_extrap();
  
  ingredients_t sub=chir.subtract_Oa2();
  
  ingredients_t evo=sub.evolve();
  
  ingredients_t ave=evo.average_equiv_momenta();
  
  ave.plot_Z("ave");
  
  // //! list of task to plot chiral extrapolation Zq
  // vector<Z_plot_task_t> Zq_plot_tasks{
  //   {&ing.Zq,&chir.Zq,&sub.Zq,&evo.Zq,string("Zq")},
  //   {&ing.Zq_sig1,&chir.Zq_sig1,&sub.Zq_sig1,&evo.Zq_sig1,"Zq_sig1"}};
  // if(use_QED) Zq_plot_tasks.push_back(make_tuple(&ing.Zq_sig1_EM,&chir.Zq_sig1_EM,&sub.Zq_sig1_EM,&evo.Zq_sig1_EM,"Zq_sig1_EM"));
  
  
  // //! list of task to print the chiral extrapolate bilinears
  // vector<Z_plot_task_t> Zbil_tasks{{&Zbil,&Zbil_chir,&Zbil_chir_sub,&Zbil_chir_sub_evolved,string("Zbil")}};
  // if(use_QED) Zbil_tasks.push_back(make_tuple(&Zbil_QED,&Zbil_QED_chir,&Zbil_QED_chir_sub,&Zbil_QED_chir_sub_evolved,"Zbil_EM"));
  
  // for(auto &p : Zbil_tasks)
  //   {
  //     //decript tuple
  //     const djvec_t &Z=(*get<0>(p));
  //     const djvec_t &Z_chir=(*get<1>(p));
  //     const djvec_t &Z_chir_sub=(*get<2>(p));
  //     const djvec_t &Z_chir_sub_evolved=(*get<3>(p));
  //     const string &tag=get<4>(p);
      
  //     for(size_t iZbil=0;iZbil<nZbil;iZbil++)
  // 	{
  // 	  //open the file
  // 	  grace_file_t out("plots/"+tag+"_Z"+Zbil_tag[iZbil]+".xmg");
  // 	  out.new_data_set();
	  
  // 	  //write mass by mass, only half of the combos
  // 	  for(size_t im1=0;im1<nm;im1++)
  // 	    for(size_t im2=im1;im2<nm;im2++)
  // 	      {
  // 		for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
  // 		  {
  // 		    djack_t y;
  // 		    y=0.0;
  // 		    for(size_t r=0;r<nr;r++) y+=Z[im_r_im_r_iZbil_indep_imom_ind({im1,r,im2,r,iZbil,indep_imom})]/nr;
		    
  // 		    out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),y.ave_err());
  // 		  }
  // 		out.new_data_set();
  // 	      }
	  
  // 	  //write chiral extrap and subtracted
  // 	  for(auto &Ztag : vector<tuple<const djvec_t*,string>>{{&Z_chir,"chir"},{&Z_chir_sub,"sub"},{&Z_chir_sub_evolved,"evo"}})
  // 	    {
  // 	      out.set_legend(get<1>(Ztag));
  // 	      for(size_t indep_imom=0;indep_imom<equiv_imoms.size();indep_imom++)
  // 		out.write_ave_err(imoms[equiv_imoms[indep_imom].first].p(L).tilde().norm2(),(*get<0>(Ztag))[iZbil_indep_imom_ind({iZbil,indep_imom})].ave_err());
  // 	      out.new_data_set();
  // 	    }
  // 	}
  //   }
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
