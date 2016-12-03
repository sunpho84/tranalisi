#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fit.hpp>
#include <functional>

using namespace placeholders;

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz(const Tpars &f0,const Tpars &b0,const Tml &ml,const Ta &a,const Tpars &adep)
{return (f0+b0*ml)*(1.0+sqr(a)*adep);}

void cont_chir_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path)
{
  //set_printlevel(3);
  
  //search max renormalized mass
  double ml_max=0;
  for(auto &data : ext_data)
    ml_max=max(ml_max,dboot_t(data.aml/a[data.ib]/z[data.ib]).ave());
  ml_max*=1.1;
  
  //parameters to fit
  MnUserParameters pars;
  vector<boot_fit_data_t> fit_data;
  
  //set a
  size_t ipar=0;
  for(size_t ia=0;ia<a.size();ia++)
    {
      pars.Add(combine("a[%zu]",ia).c_str(),a[ia][0],a[ia].err());
      fit_data.push_back(boot_fit_data_t(//numerical data
					 [a,ia]
					 (vector<double> p,int iel)
					 {return a[ia][iel];},
					 //ansatz
					 [ipar]
					 (vector<double> p,int iel)
					 {return p[ipar];},
					 //error
					 a[ia].err()));
      ipar++;
    }
  
  //set zp
  size_t offs_z=ipar;
  for(size_t iz=0;iz<z.size();iz++)
    {
      pars.Add(combine("z[%zu]",iz).c_str(),z[iz][0],z[iz].err());
      fit_data.push_back(boot_fit_data_t(//numerical data
					 [z,iz]
					 (vector<double> p,int iel)
					 {return z[iz][iel];},
					 //ansatz
					 [ipar]
					 (vector<double> p,int iel)
					 {return p[ipar];},
					 //error
					 z[iz].err()));
      ipar++;
    }
  
  //lec
  size_t if0=ipar++;
  //pars.Add("f0",0.121,0.001);
  pars.Add("f0",0.0011,0.001);
  size_t ib0=ipar++;
  //pars.Add("b0",2.57,0.01);
  pars.Add("b0",0.02,0.001);
  // size_t ic=ipar++;
  // pars.Add("c",1,1);
  // size_t iK=ipar++;
  // pars.Add("K",1,1);
  
  size_t iadep=ipar++;
  pars.Add("adep",1.8,1);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    fit_data.push_back(boot_fit_data_t(//numerical data
  				       [ext_data,idata]
  				       (vector<double> p,int iel) //dimension 2
  				       {return ext_data[idata].y[iel]/sqr(p[ext_data[idata].ib]);},
  				       //ansatz
  				       [idata,&ext_data,offs_z,iadep,ib0,if0]
  				       (vector<double> p,int iel)
  				       {
  					 double a=p[ext_data[idata].ib];
  					 double z=p[ext_data[idata].ib+offs_z];
  					 double ml=ext_data[idata].aml/a/z;
  					 return cont_chir_ansatz(p[if0],p[ib0],ml,a,p[iadep]);
  				       },
  				       //error
  				       dboot_t(ext_data[idata].y/sqr(a[ext_data[idata].ib])).err()));
  
  //! fit
  size_t iel=0;
  boot_fit_t boot_fit(fit_data,iel);
  MnMigrad migrad(boot_fit,pars);
  
  dbvec_t fit_a(a.size()),fit_z(z.size());
  dboot_t f0,b0,adep;
  dboot_t ch2;
  for(iel=0;iel<100;iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      ch2[iel]=par_min.Fval();
      
      //get back pars
      ipar=0;
      for(size_t ia=0;ia<a.size();ia++) fit_a[ia][iel]=par_min.Vec()[ipar++];
      for(size_t iz=0;iz<z.size();iz++) fit_z[iz][iel]=par_min.Vec()[ipar++];
      f0[iel]=par_min.Vec()[if0];
      b0[iel]=par_min.Vec()[ib0];
      adep[iel]=par_min.Vec()[iadep];
    }
  
  //write ch2
  cout<<"Ch2: "<<ch2.ave_err()<<" / "<<fit_data.size()-pars.Parameters().size()<<endl;
  
  //print a and z
  for(size_t ia=0;ia<a.size();ia++)
    cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<a[ia].ave_err()<<", ratio: "<<dboot_t(a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
  for(size_t iz=0;iz<z.size();iz++)
    cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<z[iz].ave_err()<<", ratio: "<<dboot_t(z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
  
  //print parameters
  cout<<"f0: "<<f0.ave_err()<<endl;
  cout<<"B0: "<<b0.ave_err()<<endl;
  cout<<"Adep: "<<adep.ave_err()<<endl;
  
  //compute physical result
  dboot_t phys_res=cont_chir_ansatz(f0, b0, ml_phys,0.0,adep);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  
  //prepare plot
  grace_file_t test(path);
  test.set_title("Continuum and chiral limit");
  test.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
  test.set_yaxis_label("$$(M^2_{\\pi^+}-M^2_{\\pi^0})/e^2 [GeV^2]");
  test.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<a.size();ib++)
    test.write_polygon(bind(cont_chir_ansatz<dboot_t,double,dboot_t>,f0,b0,_1,a[ib],adep),0,ml_max);
  //band of the continuum limit
  test.write_polygon(bind(cont_chir_ansatz<dboot_t,double,double>,f0,b0,_1,0.0,adep),0,ml_max);
  //data
  for(size_t ib=0;ib<a.size();ib++)
    {
      test.new_data_set();
       for(size_t idata=0;idata<ext_data.size();idata++)
	 if(ext_data[idata].ib==ib)
	   test<<dboot_t(ext_data[idata].aml/z[ib]/a[ib]).ave()<<" "<<dboot_t(ext_data[idata].y/sqr(a[ib])).ave_err()<<endl;
    }
  //data of the continuum-chiral limit
  test.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
  
  
}
