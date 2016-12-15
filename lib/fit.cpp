#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//linear fit

//! ansatz fit
// template <class Tpars,class Tml,class Ta>
// Tpars cont_chir_ansatz_pol(const Tpars &A,const Tpars &B,const Tpars &C,const Tml &ml,const Ta &a,const Tpars &adep)
// {
//   return A+B*ml+C*sqr(ml)+adep*sqr(a);
// }

// void cont_chir_fit_pol(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t_pol> &ext_data,const dboot_t &ml_phys,const string &path,dboot_t &output)
// {
//   //set_printlevel(3);
  
//   //search max renormalized mass
//   double ml_max=0;
//   for(auto &data : ext_data)
//     ml_max=max(ml_max,dboot_t(data.aml/a[data.ib]/z[data.ib]).ave());
//   ml_max*=1.1;
  
//   //parameters to fit
//   MnUserParameters pars;
//   vector<boot_fit_data_t> fit_data;

//   //set a
//   size_t nbeta=a.size();
//   vector<double> ipara(nbeta),iparz(nbeta);
//   for(size_t ibeta=0;ibeta<nbeta;ibeta++)
//     {
//       ipara[ibeta]=add_self_fitted_point(pars,combine("a[%zu]",ibeta),fit_data,a[ibeta]);
//       iparz[ibeta]=add_self_fitted_point(pars,combine("z[%zu]",ibeta),fit_data,z[ibeta]);
//     }
  
//   //parameters
//   size_t iA=add_fit_par(pars,"A",2,10);
//   size_t iB=add_fit_par(pars,"B",1,10);
//   size_t iC=add_fit_par(pars,"C",0,5);
//   size_t iadep=add_fit_par(pars,"adep",0.6,1);
    
//   //set data
//   for(size_t idata=0;idata<ext_data.size();idata++)
//     fit_data.push_back(boot_fit_data_t(//numerical data
//   				       [ext_data,idata]
//   				       (vector<double> p,int iel) //dimension 2
//   				       {return ext_data[idata].y[iel];},
//   				       //ansatz
//   				       [idata,&ext_data,iA,iB,iC,iadep,ipara,iparz]
//   				       (vector<double> p,int iel)
//   				       {
// 					 size_t ib=ext_data[idata].ib;
//   					 double a=p[ipara[ib]];
//   					 double z=p[iparz[ib]];
//   					 double ml=ext_data[idata].aml/a/z;
//   					 return cont_chir_ansatz_pol(p[iA],p[iB],p[iC],ml,a,p[iadep]);
//   				       },
//   				       //error
//   				       (ext_data[idata].y).err()));

//   //! fit
//   size_t iel=0;
//   boot_fit_t boot_fit(fit_data,iel);
//   MnMigrad migrad(boot_fit,pars);
  
//   dbvec_t fit_a(nbeta),fit_z(nbeta);
//   dboot_t A,B,C,adep;
//   dboot_t ch2;
//   for(iel=0;iel<ml_phys.size();iel++)
//     {
//       //minimize and print the result
//       FunctionMinimum min=migrad();
//       MinimumParameters par_min=min.Parameters();
//       ch2[iel]=par_min.Fval();
      
//       if(!min.IsValid()) CRASH("Minimization failed");
      
//       //get back pars
//       A[iel]=migrad.Value(iA);
//       B[iel]=migrad.Value(iB);
//       C[iel]=migrad.Value(iC);
//       adep[iel]=migrad.Value(iadep);
      
//       for(size_t ibeta=0;ibeta<a.size();ibeta++) fit_a[ibeta][iel]=migrad.Value(ipara[ibeta]);
//       for(size_t ibeta=0;ibeta<z.size();ibeta++) fit_z[ibeta][iel]=migrad.Value(iparz[ibeta]);
//     }
  
//   //write ch2
//   cout<<"Ch2: "<<ch2.ave_err()<<" / "<<fit_data.size()-pars.Parameters().size()<<endl;
  
//   //print a and z
//   for(size_t ia=0;ia<a.size();ia++)
//     cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<a[ia].ave_err()<<", ratio: "<<dboot_t(a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
//   for(size_t iz=0;iz<z.size();iz++)
//     cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<z[iz].ave_err()<<", ratio: "<<dboot_t(z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
  
//   //print parameters
//   cout<<"A: "<<A.ave_err()<<endl;
//   cout<<"B: "<<B.ave_err()<<endl;
//   cout<<"C: "<<C.ave_err()<<endl;
//   cout<<"Adep: "<<adep.ave_err()<<endl;
  
//   //compute physical result
//   const double a_cont=1e-5;
  
//   dboot_t phys_res=cont_chir_ansatz_pol(A,B,C,ml_phys,a_cont,adep);
//   cout<<"Physical result: "<<phys_res.ave_err()<<endl;
//   const double mk0=0.497614;
//   const double mkp=0.493677;
//   cout<<"(MK+-Mk0)^{QCD}/(-2*Deltamud): "<<(dboot_t(phys_res/(mk0+mkp))).ave_err()<<endl;
//   output=phys_res;
  
//   //prepare plot
//   grace_file_t fit_file(path);
//   fit_file.set_title("Continuum and chiral limit");
//   fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
//   fit_file.set_yaxis_label("$$(MK+-Mk0)^{QCD}/(-2*Deltamud)");
//   fit_file.set_xaxis_max(ml_max);

//   //band of the fit to individual beta
//   for(size_t ib=0;ib<a.size();ib++)
//     fit_file.write_line(bind(cont_chir_ansatz_pol<double,double,double>,A.ave(),B.ave(),C.ave(),_1,fit_a[ib].ave(),adep.ave()),1e-6,ml_max);
//   //band of the continuum limit
//   fit_file.write_polygon(bind(cont_chir_ansatz_pol<dboot_t,double,double>,A,B,C,_1,a_cont,adep),1e-6,ml_max);
//   //data
//   grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
//   for(size_t ib=0;ib<a.size();ib++)
//     {
//       fit_file.new_data_set();
//       for(size_t idata=0;idata<ext_data.size();idata++)
// 	if(ext_data[idata].ib==ib)
// 	  fit_file<<dboot_t(ext_data[idata].aml/z[ib]/a[ib]).ave()<<" "<<
// 	    (ext_data[idata].y).ave_err()<<endl;
//     }
//   //data of the continuum-chiral limit
//   fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
// }
