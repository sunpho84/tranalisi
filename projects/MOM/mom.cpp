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
  ing.plot_Z("");
  
  ingredients_t rave=ing.average_r();
  rave.plot_Z("rave");
  
  ingredients_t chir=rave.chir_extrap();
  chir.plot_Z("rave_chir");
  
  ingredients_t sub=chir.subtract_Oa2();
  sub.plot_Z("rave_chir_sub");
  
  ingredients_t evo=sub.evolve();
  evo.plot_Z("rave_chir_sub_evo");
  
  ingredients_t ave=evo.average_equiv_momenta();
  ave.plot_Z("rave_chir_sub_evo_ave");
  
  //print time statistics
  cout<<ts<<endl;
  
  return 0;
}
