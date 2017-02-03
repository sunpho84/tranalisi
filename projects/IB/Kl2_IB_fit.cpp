#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <set>

#include <fit.hpp>
#include <functional>
#include <Kl2_IB_fit.hpp>

using namespace placeholders;

namespace
{
  const double a_cont=1.0e-5;
  const double inf_vol=1.0e10;

  const double hslashc=0.197327;
  const double c_minus1_savage=-0.266596;
  const double pion_charge_radius=0.672/hslashc;
  const double kaon_charge_radius=0.56/hslashc;
  const double pion_pol_savage=7.1e-05/pow(hslashc,3.0);
  const double kaon_pol_savage=0.0;
  const double pol_savage_pi0=-1.1e-4/pow(hslashc,3.0);
  const double pol_savage_k0=-1.0e-4/pow(hslashc,3.0);
  
  const double mpi0=0.1349766;
  const double mpip=0.13957018;
  
  const double mk0=0.497614;
  const double mkp=0.493677;
}

//! compute fitted piece of FSE
template <class TL,class Ta,class TD> TD FSE_dep_L3(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_L4(const TD &D,const Ta &a,const TL &L)
{return D/(L*a*L*a*L*a*L*a);}

template <class TL,class Ta,class TD> TD FSE_dep_savage(const TD &D1,const TD &D2, const TD &D3,const TD &M,double charge_radius,double pol_savage,const Ta &a,const TL &L)
{return e2*pow(charge_radius,2.0)/(L*a*L*a*L*a)/3.0*(D1*M+D2*4*M_PI*c_minus1_savage/(L*a))+D3*8*pow(M_PI,2.0)*M*c_minus1_savage*pol_savage/(L*a*L*a*L*a*L*a);}

//! holds index and out pars
class cont_chir_fit_pars_t
{
public:
  dbvec_t ori_a,ori_z;
  dboot_t ori_f0,ori_B0;
  size_t nbeta;
  vector<size_t> ipara,iparz;
  size_t if0,iB0;
  size_t iadep,iadep_ml,iL3dep,iL4dep,iML4dep;
  size_t iC,iKPi,iKK,iK2Pi,iK2K;
  bool fitting_a,fitting_z;
  dbvec_t fit_a,fit_z;
  dboot_t fit_f0,fit_B0;
  dboot_t adep,adep_ml,L3dep,L4dep,ML4dep;
  dboot_t C,KPi,KK,K2Pi,K2K;
  
  cont_chir_fit_pars_t(size_t nbeta)
    :
    nbeta(nbeta),
    ipara(nbeta),
    iparz(nbeta),
    fitting_a(true),
    fitting_z(true),
    fit_a(nbeta),
    fit_z(nbeta)
  {}
  
  //! return a
  double get_a(const vector<double> &p,size_t ib,size_t iel) const
  {
    if(fitting_a) return p[ipara[ib]];
    else          return ori_a[ib][iel];
  }
  
  //! return z
  double get_z(const vector<double> &p,size_t ib,size_t iel) const
  {
    if(fitting_z) return p[iparz[ib]];
    else          return ori_z[ib][iel];
  }
  
  //! add a and az to be self-fitted
  void add_az_pars(const dbvec_t &a,const dbvec_t &z,boot_fit_t &boot_fit)
  {
    ori_a=a;
    ori_z=z;
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      {
	ipara[ibeta]=boot_fit.add_self_fitted_point(fit_a[ibeta],combine("a[%zu]",ibeta),a[ibeta],DO_NOT_CORRELATE);
	iparz[ibeta]=boot_fit.add_self_fitted_point(fit_z[ibeta],combine("z[%zu]",ibeta),z[ibeta],DO_NOT_CORRELATE);
      }
  }
  
  //! add adep and adep_ml
  void add_adep_pars(const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,boot_fit_t &boot_fit)
  {
    // adep, adep_ml
    iadep=boot_fit.add_fit_par(adep,"adep",adep_guess.ave,adep_guess.err);
    iadep_ml=boot_fit.add_fit_par(adep_ml,"adep_ml",adep_ml_guess.ave,adep_ml_guess.err);
  }
  
  //! add fsedep
  void add_fsedep_pars(const ave_err_t &L3dep_guess,const ave_err_t &L4dep_guess,const ave_err_t &ML4dep_guess,boot_fit_t &boot_fit)
  {
    iL3dep=boot_fit.add_fit_par(L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
    iL4dep=boot_fit.add_fit_par(L4dep,"L4dep",L4dep_guess.ave,L4dep_guess.err);
    iML4dep=boot_fit.add_fit_par(ML4dep,"ML4dep",ML4dep_guess.ave,ML4dep_guess.err);
  }
  
  //! add all common pars
  void add_common_pars(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
		       const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &L3dep_guess,const ave_err_t &L4dep_guess,const ave_err_t &ML4dep_guess,boot_fit_t &boot_fit)
  {
    add_az_pars(a,z,boot_fit);
    add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
    //f0 and B0
    ori_f0=f0;
    ori_B0=B0;
    if0=boot_fit.add_self_fitted_point(fit_f0,"f0",f0,DO_NOT_CORRELATE);
    iB0=boot_fit.add_self_fitted_point(fit_B0,"B0",B0,DO_NOT_CORRELATE);
    
    add_fsedep_pars(L3dep_guess,L4dep_guess,ML4dep_guess,boot_fit);
  }
  
  //! add the low energy constants
  void add_LEC_pars(const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &KK_guess,const ave_err_t &K2Pi_guess,const ave_err_t &K2K_guess,boot_fit_t &boot_fit)
  {
    iC=boot_fit.add_fit_par(C,"C",C_guess.ave,C_guess.err);
    iKPi=boot_fit.add_fit_par(KPi,"KPi",KPi_guess.ave,KPi_guess.err);
    iKK=boot_fit.add_fit_par(KK,"KK",KK_guess.ave,KK_guess.err);
    iK2Pi=boot_fit.add_fit_par(K2Pi,"K2Pi",K2Pi_guess.ave,K2Pi_guess.err);
    iK2K=boot_fit.add_fit_par(K2K,"K2K",K2K_guess.ave,K2K_guess.err);
  }
  
  //////////////////////////////////////////////////////////////////////
  
  void print_common_pars() const
  {
    for(size_t ia=0;ia<ori_a.size();ia++)
      cout<<"a["<<ia<<"]: "<<fit_a[ia].ave_err()<<", orig: "<<ori_a[ia].ave_err()<<", ratio: "<<dboot_t(ori_a[ia]/fit_a[ia]-1.0).ave_err()<<endl;
    for(size_t iz=0;iz<ori_z.size();iz++)
      cout<<"Zp["<<iz<<"]: "<<fit_z[iz].ave_err()<<", orig: "<<ori_z[iz].ave_err()<<", ratio: "<<dboot_t(ori_z[iz]/fit_z[iz]-1.0).ave_err()<<endl;
    //
    cout<<"f0: "<<fit_f0.ave_err()<<", orig: "<<ori_f0.ave_err()<<", ratio: "<<dboot_t(fit_f0/ori_f0-1).ave_err()<<endl;
    cout<<"B0: "<<fit_B0.ave_err()<<", orig: "<<ori_B0.ave_err()<<", ratio: "<<dboot_t(fit_B0/ori_B0-1).ave_err()<<endl;
    //
    cout<<"Adep: "<<adep.ave_err()<<endl;
    cout<<"Adep_ml: "<<adep_ml.ave_err()<<endl;
    cout<<"L3dep: "<<L3dep.ave_err()<<endl;
    cout<<"L4dep: "<<L4dep.ave_err()<<endl;
    cout<<"ML4dep: "<<ML4dep.ave_err()<<endl;
  }
  
  //! print the value
  void print_LEC_pars() const
  {
    cout<<"C: "<<C.ave_err()<<endl;
    cout<<"KPi: "<<KPi.ave_err()<<endl;
    cout<<"KK: "<<KK.ave_err()<<endl;
    cout<<"K2Pi: "<<K2Pi.ave_err()<<endl;
    cout<<"K2K: "<<K2K.ave_err()<<endl;
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
	  //make the list of volumes
	  set<size_t> L_list;
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      L_list.insert(ext_data[idata].L);
	  
	  //loop over the list of volumes
	  for(auto &L : L_list)
	    {
	      fit_file.new_data_set();
	      //put data without fse to brown
	      if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	      
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		  fit_file<<dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave()<<" "<<
		    fun_data(idata,without_with_fse,ib).ave_err()<<endl;
	    }
	}
	  
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
}

//! perform the fit to the continuum limit
void cont_chir_fit_minimize
(const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,boot_fit_t &boot_fit,double apow,double zpow,
 const function<double(const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)> &cont_chir_ansatz,bool cov_flag)
{
  //set_printlevel(3);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    boot_fit.add_point(//numerical data
		       [&ext_data,&pars,idata,apow,zpow]
		       (vector<double> p,int iel) //dimension 2
		       {return ext_data[idata].wfse[iel]*pow(pars.get_z(p,ext_data[idata].ib,iel),zpow)/pow(pars.get_a(p,ext_data[idata].ib,iel),apow);},
		       //ansatz
		       [idata,&pars,&ext_data,&cont_chir_ansatz]
		       (vector<double> p,int iel)
		       {
			 size_t ib=ext_data[idata].ib;
			 double ac=pars.get_a(p,ib,iel);
			 double zc=pars.get_z(p,ib,iel);
			 double ml=ext_data[idata].aml/ac/zc;
			 double ms=ext_data[idata].ams/ac/zc;
			 double L=ext_data[idata].L;
			 return cont_chir_ansatz(p,pars,ml,ms,ac,L);
		       },
		       //for covariance/error
		       dboot_t(ext_data[idata].wfse/pow(pars.ori_a[ext_data[idata].ib],apow)),1/*correlate*/);
  
  //! fit
  boot_fit.fit(cov_flag);
  
  //print parameters
  pars.print_common_pars();
  pars.print_LEC_pars();
}

//////////////////////////////////////////////////// pion /////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tml,class Ta>
Tpars cont_chir_ansatz_dM2Pi(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &K,const Tpars &quad_dep,const Tml &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &L4dep,const Tpars &ML4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3+16*Cf04)*M2/den*log(M2/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep*pow(M2/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_savage(L3dep,L4dep,ML4dep,sqrtM2,pion_charge_radius,pion_pol_savage,a,L);
  else                 fitted_FSE=0;
  
  return e2*sqr(f0)*(4*Cf04+chir_dep+K*M2/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2Pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,
			    const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{0.6,1.0},{0.0,10.0},{0.0,1.0},{0.0,1.0},{0.0,1.0},boot_fit);
  pars.add_LEC_pars({1.0e-5,1.0e-5},{6.0,1.0},{0.0,0.0},{0.0,10.0},{0.0,0.0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2Pi(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iL4dep],p[pars.iML4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2Pi(pars.fit_f0,B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag);
  cout<<"Physical result: "<<phys_res.ave_err()<<endl;
  cout<<"MP+-MP0: "<<(dboot_t(phys_res/(mpi0+mpip))*1000).ave_err()<<", exp: 4.5936"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2Pi<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.L4dep.ave(),pars.ML4dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_dM2Pi<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib])-
				without_with_fse*FSE_dep_savage(pars.L3dep,pars.L4dep,pars.ML4dep,dboot_t(pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)),pion_charge_radius,pion_pol_savage,pars.fit_a[ib],ext_data[idata].L));
		},
		ml_phys,phys_res,"$$(M^2_{\\pi^+}-M^2_{\\pi^0}) [GeV^2]");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_gamma ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep_pi,const Tpars &quad_dep_k,
			       const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &L4dep,const Tpars &ML4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    sqrtM2K=pow(M2K,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=M2Pi/den*log(M2Pi/pow(mu_MS,2.0))-M2K/den*log(M2K/pow(mu_MS,2.0));
  else                 chir_dep=quad_dep_pi*pow(M2Pi/den,2.0)-quad_dep_k*pow(M2K/den,2.0);
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_savage(L3dep,L4dep,ML4dep,sqrtM2K,kaon_charge_radius,kaon_pol_savage,a,L);
  else                 fitted_FSE=0;
  
  return (2+3/(4*Cf04))*(chir_dep+Kpi*M2Pi/den-Kk*M2K/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,{0.0,1.0},{0.0,0.0},{-78.0,14.0},{0.0,10.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{1.5,0.18},{-1.2,0.16},{0.0,10.0},{0.0,10.0},boot_fit);
  boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],p[pars.iK2K],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iL4dep],p[pars.iML4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon(f0,B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag);
  cout<<"epsilon_gamma: "<<phys_res.ave_err()<<", exp: 0.7"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.L4dep.ave(),pars.ML4dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-FSE_dep_savage(pars.L3dep,pars.L4dep,pars.ML4dep,dboot_t(pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)),kaon_charge_radius,kaon_pol_savage,pars.fit_a[ib],ext_data[idata].L)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\gamma");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// kaon_QED ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_dM2K_QED(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &Kk,const Tpars &quad_dep_pi,const Tpars &quad_dep_k,
				const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &L4dep,const Tpars &ML4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    sqrtM2K=pow(M2K,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=-(3+8*Cf04)*M2K/den*log(M2K/pow(mu_MS,2.0))-8*Cf04*M2Pi/den*log(M2Pi/pow(mu_MS,2.0));
  else                 chir_dep=-(quad_dep_pi*pow(M2Pi/den,2.0)+quad_dep_k*pow(M2K/den,2.0));

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_savage(L3dep,L4dep,ML4dep,sqrtM2K,kaon_charge_radius,kaon_pol_savage,a,L);
  else                 fitted_FSE=0;
  
  return e2*sqr(f0)*(4*Cf04+chir_dep+Kpi*M2Pi/den+Kk*M2K/den+disc_eff)+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_dM2K_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{2.24,1.0},{25.0,40.0},{0.23,0.05},{0.0,10.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({1.0e-4,1.0e-5},{20.0,4.0},{-50.0,10.0},{0.0,10.0},{0.0,10.0},boot_fit);
  boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,2.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_dM2K_QED(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],p[pars.iK2K],ml,ms,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iL4dep],p[pars.iML4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_dM2K_QED(f0,B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag);
  cout<<"(MK+-MK0)^{QED}: "<<(dboot_t(phys_res/(mk0+mkp))*1000).ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_dM2K_QED<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.L4dep.ave(),pars.ML4dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_dM2K_QED<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L3dep,pars.L4dep,pars.ML4dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)/sqr(pars.fit_a[ib])-
				without_with_fse*FSE_dep_savage(pars.L3dep,pars.L4dep,pars.ML4dep,dboot_t(pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)),kaon_charge_radius,kaon_pol_savage,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,"$$(M^2_{K^+}-M^2_{K^0})^{QED} [GeV^2]");
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic linear fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_linear_ansatz(const Tpars &C0,const Tpars &C1,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml;
  else                 chir_dep=0;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=sqr(a)*(adep+adep_ml*ml);
  else                 disc_eff=adep_ml*ml*sqr(a);
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_linear_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag)
{
  //set_printlevel(3);
 
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0.0,1.0},{0.0,1.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},boot_fit);
  pars.add_adep_pars({0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_fsedep_pars({0.0,1.0},{0.0,0.0},{0.0,0.0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2Pi);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_linear_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_linear_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag,with_without_FSE]
		(double x,size_t ib)
		{return cont_chir_linear_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag,with_without_FSE);},
		bind(cont_chir_linear_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,yaxis_label);
  
  return phys_res;
}

/////////////////////////////////////////////////////////////// generic quadratic fit ////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_quad_ansatz(const Tpars &C0,const Tpars &C1,const Tpars &C2,const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L3dep,const Tpars &adep_ml,const size_t an_flag,bool with_without_FSE)
{
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=C1*ml+C2*sqr(ml);
  else                 chir_dep=C1*ml;
  
  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=sqr(a)*(adep+adep_ml*ml);
  else                 disc_eff=adep_ml*ml*sqr(a);
  
  Tpars fitted_FSE;
  if(with_without_FSE==1)
    {
      if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L3(L3dep,a,L);
      else                 fitted_FSE=0;
    }
  
  return C0+chir_dep+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_quad_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,const string &yaxis_label,double apow,double zpow,size_t an_flag,bool with_without_FSE,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_az_pars(a,z,boot_fit);
  pars.add_LEC_pars({0.0,1.0},{0.0,1.0},{0.0,0.0},{0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_adep_pars({0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_fsedep_pars({0.0,1.0},{0.0,0.0},{0.0,0.0},boot_fit);
  boot_fit.fix_par(pars.iKK);
  boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,apow,zpow,[an_flag,with_without_FSE](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_quad_ansatz(p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],L,p[pars.iL3dep],p[pars.iadep_ml],an_flag,with_without_FSE);},cov_flag);
  
  dboot_t phys_res=cont_chir_quad_ansatz(pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE);
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag,with_without_FSE]
		(double x,size_t ib)
		{return cont_chir_quad_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L3dep.ave(),pars.adep_ml.ave(),an_flag,with_without_FSE);},
		bind(cont_chir_quad_ansatz<dboot_t,double,double>,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,inf_vol,pars.L3dep,pars.adep_ml,an_flag,with_without_FSE),
		[&ext_data,apow,zpow,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t((without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse)*pow(pars.fit_z[ib],zpow)/pow(pars.fit_a[ib],apow)-
				without_with_fse*FSE_dep_L3(pars.L3dep,pars.fit_a[ib],ext_data[idata].L));},
		ml_phys,phys_res,yaxis_label);
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_Pi0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_Pi0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi0,const Tpars &K,const Tpars &cubic_dep,
				   const Tm &ml,const Ta &a,const Tpars &adep,double L,const Tpars &L4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2=2*B0*ml,
    sqrtM2=pow(M2,0.5),
    den=sqr(Tpars(4*M_PI*f0));
  
  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=(3+16*Cf04)*pow(M2/(4*Cf04*den),2.0)*log(M2/pow(mu_MS,2.0))+K*pow(M2/den,2.0)/(4*Cf04)+cubic_dep*pow(M2/den,3.0);
  else                 chir_dep=K*pow(M2/den,2.0)/(4*Cf04)+cubic_dep*pow(M2/den,3.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=adep*sqr(a);

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L4(L4dep,a,L)*pol_savage_pi0*sqrtM2;
  else                 fitted_FSE=0;

  return Kpi0*(M2/(4*Cf04*den)+chir_dep)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_Pi0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{0.0,0.0},{-15.0,1.0},{0.0,0.0},{0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{10.0,1.0},{0.0,0.0},{10.0,1.0},{-25.0,1.0},boot_fit);
  boot_fit.fix_par(pars.iadep);
  //boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKK);
  //boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL3dep);
  //boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon_Pi0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],p[pars.iK2K],ml,ac,p[pars.iadep],L,p[pars.iL4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_Pi0(f0,B0,pars.C,pars.KPi,pars.K2Pi,pars.K2K,ml_phys,a_cont,pars.adep,inf_vol,pars.L4dep,pars.adep_ml,an_flag);
  cout<<"epsilon_Pi0: "<<phys_res.ave_err()<<", exp: 0.07"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_Pi0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L4dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon_Pi0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,pars.K2K,_1,a_cont,pars.adep,
		     inf_vol,pars.L4dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-pol_savage_pi0*FSE_dep_L4(pars.L4dep,pars.fit_a[ib],ext_data[idata].L)*pow(2*pars.fit_B0*ext_data[idata].aml/pars.fit_a[ib]/pars.fit_z[ib],0.5)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\Pi0");
  
  return phys_res;
}

///////////////////////////////////////////////////////// epsilon_K0 ////////////////////////////////////////////////////////////////////////////

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_epsilon_K0(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &C0,const Tpars &Kk0,const Tpars &K,const Tpars &cubic_dep,
				  const Tm &ml,const Tm &ms,const Ta &a,const Tpars &adep,double L,const Tpars &L4dep,const Tpars &adep_ml,const size_t an_flag)
{
  Tpars
    Cf04=C/sqr(sqr(f0)),
    M2Pi=2*B0*ml,
    M2K=B0*(ml+ms),
    sqrtM2K=pow(M2K,0.5),
    den=sqr(Tpars(4*M_PI*f0));

  Tpars chir_dep;
  if(chir_an(an_flag)) chir_dep=K*M2Pi/den+cubic_dep*pow(M2Pi/den,2.0)+C0*(3+16*Cf04)*M2Pi/(4*Cf04*den)*log(M2Pi/pow(mu_MS,2.0));
  else                 chir_dep=K*M2Pi/den+cubic_dep*pow(M2Pi/den,2.0);

  Tpars disc_eff;
  if(cont_an(an_flag)) disc_eff=adep*sqr(a)+adep_ml*sqr(a)*ml;
  else                 disc_eff=0;

  Tpars fitted_FSE;
  if(FSE_an(an_flag))  fitted_FSE=FSE_dep_L4(L4dep,a,L)*pol_savage_k0*sqrtM2K;
  else                 fitted_FSE=0;

  return Kk0*M2K/(4*Cf04*den)*(1+chir_dep)+disc_eff+fitted_FSE;
}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit_epsilon_K0(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,size_t an_flag,bool cov_flag)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  pars.add_common_pars(a,z,f0,B0,{-1.0,1.0},{0.0,0.0},{0.0,0.0},{0.0,1.0},{0.0,0.0},boot_fit);
  pars.add_LEC_pars({5.0e-5,1.0e-5},{1.0,0.0},{5.0,1.0},{0.0,10.0},{0.0,2.0},boot_fit);
  //boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iadep_ml);
  boot_fit.fix_par(pars.iKPi);
  //boot_fit.fix_par(pars.iK2Pi);
  //boot_fit.fix_par(pars.iK2K);
  boot_fit.fix_par(pars.iL3dep);
  //boot_fit.fix_par(pars.iL4dep);
  boot_fit.fix_par(pars.iML4dep);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_epsilon_K0(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iKK],p[pars.iK2Pi],p[pars.iK2K],ml,ms,ac,p[pars.iadep],L,p[pars.iL4dep],p[pars.iadep_ml],an_flag);},cov_flag);
  
  dboot_t phys_res=cont_chir_ansatz_epsilon_K0(f0,B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,ml_phys,ms_phys,a_cont,pars.adep,inf_vol,pars.L4dep,pars.adep_ml,an_flag);
  cout<<"epsilon_K0: "<<phys_res.ave_err()<<", exp: 0.3"<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,&ms_phys,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_epsilon_K0<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.KK.ave(),pars.K2Pi.ave(),pars.K2K.ave(),x,ms_phys.ave(),
		     pars.fit_a[ib].ave(),pars.adep.ave(),inf_vol,pars.L4dep.ave(),pars.adep_ml.ave(),an_flag);},
		bind(cont_chir_ansatz_epsilon_K0<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.KK,pars.K2Pi,pars.K2K,_1,ms_phys.ave(),a_cont,pars.adep,
		     inf_vol,pars.L4dep,pars.adep_ml,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?(ext_data[idata].wfse-pol_savage_k0*FSE_dep_L4(pars.L4dep,pars.fit_a[ib],ext_data[idata].L)*pow(pars.fit_B0*(ext_data[idata].aml+ext_data[idata].ams)/pars.fit_a[ib]/pars.fit_z[ib],0.5)):ext_data[idata].wofse);},
		ml_phys,phys_res,"$$\\varepsilon_\\K0");
  
  return phys_res;
}

