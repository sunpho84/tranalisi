#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fit.hpp>
#include <functional>
#include <Kl2_IB_fit.hpp>

using namespace placeholders;

namespace
{
  double a_cont=1e-5;
  double inf_vol=1e10;
  
  double mpi0=0.1349766;
  double mpip=0.13957018;
}

//! compute fitted piece of FSE
template <class TL,class Ta,class TD> TD FSE_dep(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a);}

//! holds index and out pars
class cont_chir_fit_pars_t
{
  size_t nbeta;
public:
  vector<size_t> ipara,iparz;
  size_t if0,iB0;
  size_t iadep,iadep_ml,iL3dep;
  size_t iC,iKPi,iKK;
  dbvec_t fit_a,fit_z;
  dboot_t fit_f0,fit_B0;
  dboot_t adep,adep_ml,L3dep;
  dboot_t C,KPi,KK;
  
  cont_chir_fit_pars_t(size_t nbeta) : nbeta(nbeta),ipara(nbeta),iparz(nbeta),fit_a(nbeta),fit_z(nbeta) {}
  
  //! add all common pars
  void add_common_pars(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
		       const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &L3dep_guess,boot_fit_t &boot_fit)
  {
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      {
	ipara[ibeta]=boot_fit.add_self_fitted_point(fit_a[ibeta],combine("a[%zu]",ibeta),a[ibeta]);
	iparz[ibeta]=boot_fit.add_self_fitted_point(fit_z[ibeta],combine("z[%zu]",ibeta),z[ibeta]);
      }
    //f0 and B0
    if0=boot_fit.add_self_fitted_point(fit_f0,"f0",f0);
    iB0=boot_fit.add_self_fitted_point(fit_B0,"B0",B0);
    // adep, adep_ml and L3 dep
    iadep=boot_fit.add_fit_par(adep,"adep",adep_guess.ave,adep_guess.err);
    iadep_ml=boot_fit.add_fit_par(adep_ml,"adep_ml",adep_ml_guess.ave,adep_ml_guess.err);
    iL3dep=boot_fit.add_fit_par(L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  }
  
  //! add the low energy constants
  void add_LEC_pars(const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &KK_guess,boot_fit_t &boot_fit)
  {
    iC=boot_fit.add_fit_par(C,"C",C_guess.ave,C_guess.err);
    iKPi=boot_fit.add_fit_par(KPi,"KPi",KPi_guess.ave,KPi_guess.err);
    iKK=boot_fit.add_fit_par(KK,"KK",KK_guess.ave,KK_guess.err);
  }
  
  //////////////////////////////////////////////////////////////////////
  
  void print_common_pars(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0)
  {
    for(size_t ia=0;ia<a.size();ia++)
      cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<a[ia].ave_err()<<", ratio: "<<dboot_t(a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
    for(size_t iz=0;iz<z.size();iz++)
      cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<z[iz].ave_err()<<", ratio: "<<dboot_t(z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
    //
    cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/f0-1).ave_err()<<endl;
    cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/B0-1).ave_err()<<endl;
    //
    cout<<"Adep: "<<adep.ave_err()<<endl;
    cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
    cout<<"L3dep: "<<L3dep.ave_err()<<endl;
  }
  
  //! print the value
  void print_LEC_pars()
  {
    cout<<"C: "<<C.ave_err()<<endl;
    cout<<"KPi: "<<KPi.ave_err()<<endl;
    cout<<"KK: "<<KK.ave_err()<<endl;
  }
};

//! plot the continuum-chiral extrapolation
void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label)
{
  //search max renormalized mass
  double ml_max=0;
  for(auto &data : ext_data)
    ml_max=max(ml_max,dboot_t(data.aml/pars.fit_a[data.ib]/pars.fit_z[data.ib]).ave());
  ml_max*=1.1;
  
  //prepare plot
  grace_file_t fit_file(path);
  fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<pars.fit_a.size();ib++)
	{
	  fit_file.new_data_set();
	  //put data without fse to brown
	  if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	  
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      fit_file<<dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave()<<" "<<
		fun_data(idata,without_with_fse,ib).ave_err()<<endl;
	}
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
}

//////////////////////////////////////////////////// pion /////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tml &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  if(chir_an) chir_log=-(3+16*Cf04)*M2/den*log(M2/sqr(mu));
  else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return e2*sqr(f0)*(4*Cf04+chir_log+K*M2/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
cont_chir_fit_pars_t cont_chir_fit_dM2Pi_minimize(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,bool chir_an)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,{0.6,1.0},{0.0,10.0},{0.0,1.0},boot_fit);
  pars.add_LEC_pars({1e-5,1e-5},{6,1},{0,0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    boot_fit.add_point(//numerical data
		       [&ext_data,&pars,idata]
		       (vector<double> p,int iel) //dimension 2
		       {return ext_data[idata].wfse[iel]/sqr(p[pars.ipara[ext_data[idata].ib]]);},
		       //ansatz
		       [idata,&pars,&ext_data,chir_an]
		       (vector<double> p,int iel)
		       {
			 size_t ib=ext_data[idata].ib;
			 double a=p[pars.ipara[ib]];
			 double z=p[pars.iparz[ib]];
			 double ml=ext_data[idata].aml/a/z;
			 double L=ext_data[idata].L;
			 return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],ml,a,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],chir_an);
		       },
		       //error
		       dboot_t(ext_data[idata].wfse/sqr(a[ext_data[idata].ib])).err());
  
  //! fit
  boot_fit.fit();
  
  //print parameters
  pars.print_common_pars(a,z,f0,B0);
  pars.print_LEC_pars();
  
  return pars;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
		      const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an)
{
  cont_chir_fit_pars_t pars=cont_chir_fit_dM2Pi_minimize(a,z,f0,B0,ext_data,chir_an);
  
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,B0,pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2Pi<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,
		     pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib])-
				without_with_fse*FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,"$$(M^2_{\\pi^+}-M^2_{\\pi^0}) [GeV^2]");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,
			       const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const bool chir_an)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_log;
  // if(chir_an)
  chir_log=M2Pi/den*log(M2Pi/sqr(mu))-M2K/den*log(M2K/sqr(mu));
  // else        chir_log=0;
  
  Tpars fitted_FSE=FSE_dep(L3dep,a,L);
  Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
  return (2+3/(4*Cf04))*(chir_log+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
cont_chir_fit_pars_t cont_chir_fit_epsilon_minimize(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,bool chir_an)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,{0,1},{0.0,0.0},{-78,14},boot_fit);
  pars.add_LEC_pars({5e-5,1e-5},{1.5,0.18},{-1.2,0.16},boot_fit);
  boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    boot_fit.add_point(//numerical data
		       [ext_data,idata]
		       (vector<double> p,int iel) //dimension 2
		       {return ext_data[idata].wfse[iel];},
		       //ansatz
				       [idata,&pars,&ext_data,chir_an]
		       (vector<double> p,int iel)
		       {
			 size_t ib=ext_data[idata].ib;
			 double a=p[pars.ipara[ib]];
			 double z=p[pars.iparz[ib]];
			 double ml=ext_data[idata].aml/a/z;
			 double ms=ext_data[idata].ams/a/z;
			 double L=ext_data[idata].L;
			 return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],ml,ms,a,p[pars.iadep],L,
							 p[pars.iL3dep],p[pars.iadep_ml],chir_an);
		       },
		       //error
		       ext_data[idata].wfse.err());
  
  //! fit
  boot_fit.fit();
  
  //print parameters
  pars.print_common_pars(a,z,f0,B0);
  pars.print_LEC_pars();
  
  return pars;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			      const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an)
{
  cont_chir_fit_pars_t pars=cont_chir_fit_epsilon_minimize(a,z,f0,B0,ext_data,chir_an);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(f0,B0,pars.C,pars.KPi,pars.KK,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,chir_an);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,chir_an]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),chir_an);},
		bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.adep_ml,chir_an),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(-(without_with_fse-1)*ext_data[idata].wofse+
				without_with_fse*(ext_data[idata].wfse-FSE_dep(pars.L3dep,pars.fit_a[ib],ext_data[idata].L)));},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma");
  
  return phys_res;
}

// //kaon_QED

// //! ansatz fit
// template <class Tpars,class Tm,class Ta>
// Tpars cont_chir_ansatz_k(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &D,const Tpars &adep_ml,const bool chir_an)
// {
//   Tpars
//     Cf04=C/sqr(sqr(f0)),
//     M2Pi=2*B0*ml,
//     M2K=B0*(ml+ms),
//     den=sqr(Tpars(4*M_PI*f0));
  
//   Tpars chir_log;
//   if(chir_an) chir_log=-(3+8*Cf04)*M2K/den*log(M2K/sqr(mu))-8*Cf04*M2Pi/den*log(M2Pi/sqr(mu));
//   else        chir_log=0;
  
//   Tpars fitted_FSE=FSE_dep(D,a,L);
//   Tpars disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  
//   return e2*sqr(f0)*(4*Cf04+chir_log+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
// }

// void cont_chir_fit_k(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_k> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an)
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
  
//   //lec
//   size_t if0=add_self_fitted_point(pars,"f0",fit_data,f0);
//   size_t iB0=add_self_fitted_point(pars,"B0",fit_data,B0);
//   size_t iC=add_fit_par(pars,"C",1e-4,1e-5);
//   size_t iKpi=add_fit_par(pars,"Kpi",20,4);
//   size_t iKk=add_fit_par(pars,"Kk",-50,10);
  
//   size_t iadep=add_fit_par(pars,"adep",2.24,1);
//   size_t iD=add_fit_par(pars,"D",0.23,0.05);
//   size_t iadep_ml=add_fit_par(pars,"adep_ml",25,40);
  
//   //set data
//   for(size_t idata=0;idata<ext_data.size();idata++)
//     fit_data.push_back(boot_fit_data_t(//numerical data
//   				       [ext_data,idata]
//   				       (vector<double> p,int iel) //dimension 2
//   				       {return (ext_data[idata].y[iel]-ext_data[idata].fse[iel])/sqr(p[2*ext_data[idata].ib+0]);},
//   				       //ansatz
//   				       [idata,&ext_data,iadep_ml,iD,iadep,iB0,if0,iC,iKpi,iKk,chir_an,ipara,iparz]
//   				       (vector<double> p,int iel)
//   				       {
// 					 size_t ib=ext_data[idata].ib;
//   					 double a=p[ipara[ib]];
//   					 double z=p[iparz[ib]];
//   					 double ml=ext_data[idata].aml/a/z;
// 					 double ms=ext_data[idata].ams/a/z;
// 					 double L=ext_data[idata].L;
//   					 return cont_chir_ansatz_k(p[if0],p[iB0],p[iC],p[iKpi],p[iKk],ml,ms,a,p[iadep],L,p[iD],p[iadep_ml],chir_an);
//   				       },
//   				       //error
//   				       dboot_t((ext_data[idata].y-ext_data[idata].fse)/sqr(a[ext_data[idata].ib])).err()));

//   //! fit
//   size_t iel=0;
//   boot_fit_t boot_fit(fit_data,iel);
//   MnMigrad migrad(boot_fit,pars);
  
//   dbvec_t fit_a(nbeta),fit_z(nbeta);
//   dboot_t fit_f0,fit_B0,C,Kpi,Kk,adep,D,adep_ml;
//   dboot_t ch2;
//   for(iel=0;iel<ml_phys.size();iel++)
//     {
//       //minimize and print the result
//       FunctionMinimum min=migrad();
//       MinimumParameters par_min=min.Parameters();
//       ch2[iel]=par_min.Fval();
      
//       if(!min.IsValid()) CRASH("Minimization failed");
      
//       //get back pars
//       fit_f0[iel]=migrad.Value(if0);
//       fit_B0[iel]=migrad.Value(iB0);
//       C[iel]=migrad.Value(iC);
//       Kpi[iel]=migrad.Value(iKpi);
//       Kk[iel]=migrad.Value(iKk);
//       adep[iel]=migrad.Value(iadep);
//       D[iel]=migrad.Value(iD);
//       adep_ml[iel]=migrad.Value(iadep_ml);
      
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
//   cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/f0-1).ave_err()<<endl;
//   cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/B0-1).ave_err()<<endl;
//   cout<<"C: "<<C.ave_err()<<endl;
//   cout<<"Kpi: "<<Kpi.ave_err()<<endl;
//   cout<<"Kk: "<<Kk.ave_err()<<endl;
//   cout<<"Adep: "<<adep.ave_err()<<endl;
//   cout<<"D: "<<D.ave_err()<<endl;
//   cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
  
//   //compute physical result
//   const double a_cont=1e-5;
//   const double inf_vol=1e10;
  
//   dboot_t phys_res=cont_chir_ansatz_k(f0,B0,C,Kpi,Kk,ml_phys,ms_phys,a_cont,adep,inf_vol,D,adep_ml,chir_an);
//   cout<<"Physical result: "<<phys_res.ave_err()<<endl;
//   const double mk0=0.497614;
//   const double mkp=0.493677;
//   cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;
  
//   //prepare plot
//   grace_file_t fit_file(path);
//   fit_file.set_title("Continuum and chiral limit");
//   fit_file.set_xaxis_label("$$ml^{\\overline{MS},2 GeV} [GeV]");
//   fit_file.set_yaxis_label("$$(M^2_{\\K^+}-M^2_{\\K^0})^{QED} [GeV^2]");
//   fit_file.set_xaxis_max(ml_max);

//   //band of the fit to individual beta
//   for(size_t ib=0;ib<a.size();ib++)
//     fit_file.write_line(bind(cont_chir_ansatz_k<double,double,double>,fit_f0.ave(),fit_B0.ave(),C.ave(),Kpi.ave(),Kk.ave(),_1,ms_phys.ave(),fit_a[ib].ave(),adep.ave(),inf_vol,D.ave(),adep_ml.ave(),chir_an),1e-6,ml_max);
//   //band of the continuum limit
//   fit_file.write_polygon(bind(cont_chir_ansatz_k<dboot_t,double,double>,fit_f0,fit_B0,C,Kpi,Kk,_1,ms_phys.ave(),a_cont,adep,inf_vol,D,adep_ml,chir_an),1e-6,ml_max);
//   //data without and with fse
//   grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
//   for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
//     {
//       for(size_t ib=0;ib<a.size();ib++)
// 	{
// 	  fit_file.new_data_set();
// 	  //put data without fse to brown
// 	  if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	  
// 	  for(size_t idata=0;idata<ext_data.size();idata++)
// 	    if(ext_data[idata].ib==ib)
// 	      fit_file<<dboot_t(ext_data[idata].aml/z[ib]/a[ib]).ave()<<" "<<
// 		dboot_t((ext_data[idata].y-without_with_fse*ext_data[idata].fse)/sqr(fit_a[ib])-without_with_fse*FSE_dep(D,fit_a[ib],ext_data[idata].L)).ave_err()<<endl;
// 	}
//       //put back colors for data with fse
//       if(without_with_fse==0) fit_file.reset_cur_col();
//     }
//   //data of the continuum-chiral limit
//   fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
// }

