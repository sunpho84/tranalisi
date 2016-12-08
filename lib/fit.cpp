#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fit.hpp>
#include <functional>

using namespace placeholders;

//! compute fitted piece of FSE
template <class TL,class Ta,class TD> TD FSE_dep(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a);}

//pion

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tml &ml,const Ta &a,const Tpars &adep,double L,const Tpars &D,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=-(3+16*Cf04)*M2/den*log(M2/sqr(mu));
  else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(D,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return e2*sqr(f0)*(4*Cf04+chir_log+K*M2/den+disc_eff)+fitted_FSE;
}

void cont_chir_fit_pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_pi> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an)
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
  size_t nbeta=a.size();
  vector<double> ipara(nbeta),iparz(nbeta);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      ipara[ibeta]=add_self_fitted_point(pars,combine("a[%zu]",ibeta),fit_data,a[ibeta]);
      iparz[ibeta]=add_self_fitted_point(pars,combine("z[%zu]",ibeta),fit_data,z[ibeta]);
    }
  
  //lec
  size_t if0=add_self_fitted_point(pars,"f0",fit_data,f0);
  size_t iB0=add_self_fitted_point(pars,"B0",fit_data,B0);
  size_t iC=add_fit_par(pars,"C",1e-5,1e-5);
  size_t iK=add_fit_par(pars,"K",6,10);
  
  size_t iadep=add_fit_par(pars,"adep",0.6,1);
  size_t iD=add_fit_par(pars,"D",0,1);
  size_t iadep_ml=add_fit_par(pars,"adep_ml",0,10);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    fit_data.push_back(boot_fit_data_t(//numerical data
  				       [ext_data,idata]
  				       (vector<double> p,int iel) //dimension 2
  				       {return (ext_data[idata].y[iel]-ext_data[idata].fse[iel])/sqr(p[2*ext_data[idata].ib+0]);},
  				       //ansatz
  				       [idata,&ext_data,iadep_ml,iD,iadep,iB0,if0,iC,iK,chir_an,ipara,iparz]
  				       (vector<double> p,int iel)
  				       {
					 size_t ib=ext_data[idata].ib;
  					 double a=p[ipara[ib]];
  					 double z=p[iparz[ib]];
  					 double ml=ext_data[idata].aml/a/z;
					 double L=ext_data[idata].L;
  					 return cont_chir_ansatz_pi(p[if0],p[iB0],p[iC],p[iK],ml,a,p[iadep],L,p[iD],p[iadep_ml],chir_an);
  				       },
  				       //error
  				       dboot_t((ext_data[idata].y-ext_data[idata].fse)/sqr(a[ext_data[idata].ib])).err()));

  //! fit
  size_t iel=0;
  boot_fit_t boot_fit(fit_data,iel);
  MnMigrad migrad(boot_fit,pars);
  
  dbvec_t fit_a(nbeta),fit_z(nbeta);
  dboot_t fit_f0,fit_B0,C,K,adep,D,adep_ml;
  dboot_t ch2;
  for(iel=0;iel<ml_phys.size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      ch2[iel]=par_min.Fval();
      
      if(!min.HasValidParameters()) CRASH("Minimization failed");
      
      //get back pars
      fit_f0[iel]=par_min.Vec()[if0];
      fit_B0[iel]=par_min.Vec()[iB0];
      C[iel]=par_min.Vec()[iC];
      K[iel]=par_min.Vec()[iK];
      adep[iel]=par_min.Vec()[iadep];
      D[iel]=par_min.Vec()[iD];
      adep_ml[iel]=par_min.Vec()[iadep_ml];
      
      for(size_t ibeta=0;ibeta<a.size();ibeta++) fit_a[ibeta][iel]=par_min.Vec()[ipara[ibeta]];
      for(size_t ibeta=0;ibeta<z.size();ibeta++) fit_z[ibeta][iel]=par_min.Vec()[iparz[ibeta]];
    }
  
  //write ch2
  cout<<"Ch2: "<<ch2.ave_err()<<" / "<<fit_data.size()-pars.Parameters().size()<<endl;
  
  //print a and z
  for(size_t ia=0;ia<a.size();ia++)
    cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<a[ia].ave_err()<<", ratio: "<<dboot_t(a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
  for(size_t iz=0;iz<z.size();iz++)
    cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<z[iz].ave_err()<<", ratio: "<<dboot_t(z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
  
  //print parameters
  cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/f0-1).ave_err()<<endl;
  cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/B0-1).ave_err()<<endl;
  cout<<"C: "<<C.ave_err()<<endl;
  cout<<"K: "<<K.ave_err()<<endl;
  cout<<"Adep: "<<adep.ave_err()<<endl;
  cout<<"D: "<<D.ave_err()<<endl;
  cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
  
  //compute physical result
  const double a_cont=1e-5;
  const double inf_vol=1e10;
  
  dboot_t phys_res=cont_chir_ansatz_pi(f0,B0,C,K,ml_phys,a_cont,adep,inf_vol,D,adep_ml,chir_an);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  //prepare plot
  grace_file_t fit_file(path);
  fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
  fit_file.set_yaxis_label("$$(M^2_{\\pi^+}-M^2_{\\pi^0}) [GeV^2]");
  fit_file.set_xaxis_max(ml_max);

  //band of the fit to individual beta
  for(size_t ib=0;ib<a.size();ib++)
    fit_file.write_line(bind(cont_chir_ansatz_pi<double,double,double>,fit_f0.ave(),fit_B0.ave(),C.ave(),K.ave(),_1,fit_a[ib].ave(),adep.ave(),inf_vol,D.ave(),adep_ml.ave(),chir_an),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(bind(cont_chir_ansatz_pi<dboot_t,double,double>,fit_f0,fit_B0,C,K,_1,a_cont,adep,inf_vol,D,adep_ml,chir_an),1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<a.size();ib++)
	{
	  fit_file.new_data_set();
	  //put data without fse to brown
	  if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	  
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      fit_file<<dboot_t(ext_data[idata].aml/z[ib]/a[ib]).ave()<<" "<<
		dboot_t((ext_data[idata].y-without_with_fse*ext_data[idata].fse)/sqr(fit_a[ib])-without_with_fse*FSE_dep(D,fit_a[ib],ext_data[idata].L)).ave_err()<<endl;
	}
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
}
/*
//epsilon_gamma

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &D,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=M2Pi/den*log(M2Pi/sqr(mu))-M2K/den*log(M2K/sqr(mu));
  else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(D,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return (2+3/(4*Cf04))*(chir_log+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
}

void cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_epsilon> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an)
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
  size_t nbeta=a.size();
  vector<double> ipara(nbeta),iparz(nbeta);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      ipara[ibeta]=add_self_fitted_point(pars,combine("a[%zu]",ibeta),fit_data,a[ibeta]);
      iparz[ibeta]=add_self_fitted_point(pars,combine("z[%zu]",ibeta),fit_data,z[ibeta]);
    }
  
  //lec
  size_t if0=add_self_fitted_point(pars,"f0",fit_data,f0);
  size_t iB0=add_self_fitted_point(pars,"B0",fit_data,B0);
  size_t iC=add_fit_par(pars,"C",1e-5,1e-5);
  size_t iKpi=add_fit_par(pars,"Kpi",0,1);
  size_t iKk=add_fit_par(pars,"Kk",0,1);
  
  size_t iadep=add_fit_par(pars,"adep",0.6,1);
  size_t iD=add_fit_par(pars,"D",0,1);
  size_t iadep_ml=add_fit_par(pars,"adep_ml",0,10);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    fit_data.push_back(boot_fit_data_t(//numerical data
  				       [ext_data,idata]
  				       (vector<double> p,int iel) //dimension 2
  				       {return ext_data[idata].wfse[iel];},
  				       //ansatz
  				       [idata,&ext_data,iadep_ml,iD,iadep,iB0,if0,iC,iKpi,iKk,chir_an,ipara,iparz]
  				       (vector<double> p,int iel)
  				       {
					 size_t ib=ext_data[idata].ib;
  					 double a=p[ipara[ib]];
  					 double z=p[iparz[ib]];
  					 double ml=ext_data[idata].aml/a/z;
					 double ms=ext_data[idata].ams/a/z;
					 double L=ext_data[idata].L;
  					 return cont_chir_ansatz_epsilon(p[if0],p[iB0],p[iC],p[iKpi],p[iKk],ml,ms,a,p[iadep],L,p[iD],p[iadep_ml],chir_an);
  				       },
  				       //error
  				       (ext_data[idata].wfse).err()));

  //! fit
  size_t iel=0;
  boot_fit_t boot_fit(fit_data,iel);
  MnMigrad migrad(boot_fit,pars);
  
  dbvec_t fit_a(nbeta),fit_z(nbeta);
  dboot_t fit_f0,fit_B0,C,Kpi,Kk,adep,D,adep_ml;
  dboot_t ch2;
  for(iel=0;iel<ml_phys.size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      ch2[iel]=par_min.Fval();
      
      if(!min.HasValidParameters()) CRASH("Minimization failed");
      
      //get back pars
      fit_f0[iel]=par_min.Vec()[if0];
      fit_B0[iel]=par_min.Vec()[iB0];
      C[iel]=par_min.Vec()[iC];
      Kpi[iel]=par_min.Vec()[iKpi];
      Kk[iel]=par_min.Vec()[iKk];
      adep[iel]=par_min.Vec()[iadep];
      D[iel]=par_min.Vec()[iD];
      adep_ml[iel]=par_min.Vec()[iadep_ml];
      
      for(size_t ibeta=0;ibeta<a.size();ibeta++) fit_a[ibeta][iel]=par_min.Vec()[ipara[ibeta]];
      for(size_t ibeta=0;ibeta<z.size();ibeta++) fit_z[ibeta][iel]=par_min.Vec()[iparz[ibeta]];
    }
  
  //write ch2
  cout<<"Ch2: "<<ch2.ave_err()<<" / "<<fit_data.size()-pars.Parameters().size()<<endl;
  
  //print a and z
  for(size_t ia=0;ia<a.size();ia++)
    cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<a[ia].ave_err()<<", ratio: "<<dboot_t(a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
  for(size_t iz=0;iz<z.size();iz++)
    cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<z[iz].ave_err()<<", ratio: "<<dboot_t(z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
  
  //print parameters
  cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/f0-1).ave_err()<<endl;
  cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/B0-1).ave_err()<<endl;
  cout<<"C: "<<C.ave_err()<<endl;
  cout<<"Kpi: "<<Kpi.ave_err()<<endl;
  cout<<"Kk: "<<Kk.ave_err()<<endl;
  cout<<"Adep: "<<adep.ave_err()<<endl;
  cout<<"D: "<<D.ave_err()<<endl;
  cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
  
  //compute physical result
  const double a_cont=1e-5;
  const double inf_vol=1e10;
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(f0,B0,C,Kpi,Kk,ml_phys,ms_phys,a_cont,adep,inf_vol,D,adep_ml,chir_an);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;
  
  //prepare plot
  grace_file_t fit_file(path);
  fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
  fit_file.set_yaxis_label("$$epsilon_gamma");
  fit_file.set_xaxis_max(ml_max);

  //band of the fit to individual beta
  for(size_t ib=0;ib<a.size();ib++)
    fit_file.write_line(bind(cont_chir_ansatz_epsilon<double,double,double>,fit_f0.ave(),fit_B0.ave(),C.ave(),Kpi.ave(),Kk.ave(),_1,ms_phys.ave(),fit_a[ib].ave(),adep.ave(),inf_vol,D.ave(),adep_ml.ave(),chir_an),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,fit_f0,fit_B0,C,Kpi,Kk,_1,ms_phys,a_cont,adep,inf_vol,D,adep_ml,chir_an),1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<a.size();ib++)
	{
	  fit_file.new_data_set();
	  //put data without fse to brown
	  if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	  
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      fit_file<<dboot_t(ext_data[idata].aml/z[ib]/a[ib]).ave()<<" "<<
		dboot_t(-(without_with_fse-1)*ext_data[idata].wofse+without_with_fse*(ext_data[idata].wfse-FSE_dep(D,fit_a[ib],ext_data[idata].L))).ave_err()<<endl;
	}
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
}
*/
