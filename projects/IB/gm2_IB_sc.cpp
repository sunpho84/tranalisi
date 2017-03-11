#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr
//const int im=istrange;
const int im=icharm;

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz_LO(const Tpars &f0,const Tpars &B0,const Tpars &C,const Tpars &Kpi,const Tpars &K2pi,const Tm &ml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L3dep,const size_t an_flag)
{
    Tpars
      M2=2*B0*ml,M=sqrt(M2),
      den=sqr(Tpars(4*M_PI*f0)),
      xi=M2/den;
    
    return C*(1+Kpi*xi+K2pi*xi*xi+a*a*(adep+ml*adep_ml)+L3dep*xi*exp(-M*L)/(M*L));
}

//! perform the fit to the continuum limit of LO
dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list,const ave_err_t &adep_guess,const ave_err_t &adep_ml_guess,const ave_err_t &C_guess,const ave_err_t &KPi_guess,const ave_err_t &L3dep_guess,const string &title)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  pars.add_common_pars(a,z,f0,B0,adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C",C_guess.ave,C_guess.err);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",KPi_guess.ave,KPi_guess.err);
  pars.iK2Pi=boot_fit.add_fit_par(pars.K2Pi,"K2Pi",KPi_guess.ave,KPi_guess.err);
  pars.iL3dep=boot_fit.add_fit_par(pars.L3dep,"L3dep",L3dep_guess.ave,L3dep_guess.err);
  
  //boot_fit.fix_par(pars.iadep);
  if(FSE_an(an_flag)) boot_fit.fix_par(pars.iL3dep);
  boot_fit.fix_par(pars.iKPi);
  if(cont_an(an_flag)) boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz_LO(p[pars.if0],p[pars.iB0],p[pars.iC],p[pars.iKPi],p[pars.iK2Pi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL3dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz_LO(pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz_LO<double,double,double>
		    (pars.fit_f0.ave(),pars.fit_B0.ave(),pars.C.ave(),pars.KPi.ave(),pars.K2Pi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L3dep.ave(),an_flag);},
		bind(cont_chir_ansatz_LO<dboot_t,double,double>,pars.fit_f0,pars.fit_B0,pars.C,pars.KPi,pars.K2Pi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L3dep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse// -without_with_fse*pars.L2dep/sqr(ext_data[idata].L)
				);},
		ml_phys,phys_res,title,beta_list);
  
  return phys_res;
}

dboot_t cont_chir_fit_LO(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  const ave_err_t adep_guess={0,0};
  const ave_err_t adep_ml_guess={0,0};
  const ave_err_t C_guess={0,0};
  const ave_err_t KPi_guess={0,0};
  const ave_err_t L3dep_guess={0,0};
  const string title="$$a_\\mu";
  
  return cont_chir_fit_LO(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,L3dep_guess,title);
}

// dboot_t cont_chir_fit_QED(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_guess={0,0};
//   const ave_err_t adep_ml_guess={0,0};
//   const ave_err_t C_guess={0,0};
//   const ave_err_t KPi_guess={0,0};
//   const ave_err_t L3dep_guess={0,0};
//   const string title="$$\\delta a_\\mu";
  
//   return cont_chir_fit(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,L3dep_guess,title);
// }

// dboot_t cont_chir_fit_ratio(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
// {
//   const ave_err_t adep_guess={0,0};
//   const ave_err_t adep_ml_guess={0,0};
//   const ave_err_t C_guess={0,0};
//   const ave_err_t KPi_guess={0,0};
//   const ave_err_t L3dep_guess={0,0};
//   const string title="$$\\delta a_\\mu/a_\\mu";
  
//   return cont_chir_fit(a,z,f0,B0,ext_data,ml_phys,path,an_flag,cov_flag,beta_list,adep_guess,adep_ml_guess,C_guess,KPi_guess,L3dep_guess,title);
// }

void perform_analysis(const dbvec_t &data,const string &name)
{
  vector<ave_err_t> ave=ave_analyses(data);
  cout<<name<<": "<<stat_analysis(ave)<<" "<<syst_analysis(ave)<<endl;
   for(size_t i=0;i<ave.size();i++)
     cout<<"  "<<ave[i]<<endl;
}

int main(int narg,char **arg)
{
  gm2_initialize(narg,arg);
  
  djvec_t bare_LO(ens_data.size());
  djvec_t bare_QED(ens_data.size());
  djvec_t ratio(ens_data.size());
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t &ens=ens_data[iens];
      size_t TH=ens.T/2;
      cout<<"----------------------- "<<iens<<" "<<ens.path<<" ---------------------"<<endl;
      
      djack_t deltam_cr=compute_deltam_cr(ens,ilight);
      
      //load LO for PP
      djvec_t PP_LO=read_PP("00",ens,im,1,RE);
      djack_t Z_P,M_P;
      two_pts_fit(Z_P,M_P,PP_LO,TH,ens.tmin[im], ens.tmax[im],combine("%s/plots/PP_LO_emass.xmg",ens.path.c_str()));
      
      //load LO and QED for VV
      djvec_t VV_LO=read_VV("00",ens,im,1,RE),eff_VV_LO=effective_mass(VV_LO,TH);
      djvec_t VV_QED=read_QED_VV(ens,1,im,deltam_cr,VV_LO),VV_rat=VV_QED/VV_LO;
      djack_t Z_V,M_V,A_V,SL_V;
      two_pts_with_ins_ratio_fit(Z_V,M_V,A_V,SL_V,VV_LO,VV_QED,TH,ens.tmin[im],ens.tmax[im],
				 combine("%s/plots/VV_LO.xmg",ens.path.c_str()),combine("%s/plots/VV_QED.xmg",ens.path.c_str()));
      
      //lattice spacing obtained from V
      djack_t resc_a=M_V/M_V_phys[im];
      
      size_t upto=TH-DT;
      //cout<<"Upto: "<<upto;
      djack_t bare_LO_correl=integrate_corr_times_kern_up_to(VV_LO,ens.T,resc_a,im,upto);
      //cout<<", after: "<<upto<<endl;
      djack_t bare_LO_remaind=integrate_LO_reco_from(Z_V,M_V,resc_a,im,upto);
      bare_LO[iens]=bare_LO_correl+bare_LO_remaind;
      compare_LO_num_reco(combine("%s/plots/kern_LO_num_reco.xmg",ens.path.c_str()),VV_LO,Z_V,M_V,resc_a);
      cout<<"bare amu: "<< bare_LO_correl.ave_err()<<" + "<<bare_LO_remaind.ave_err()<<" = "<<bare_LO[iens].ave_err()<<endl;
      
      djack_t bare_QED_correl=integrate_corr_times_kern_up_to(VV_QED,ens.T,resc_a,im,upto);
      djack_t bare_QED_remaind=integrate_QED_reco_from(A_V,Z_V,SL_V,M_V,resc_a,im,upto);
      bare_QED[iens]=bare_QED_correl+bare_QED_remaind;
      compare_QED_num_reco(combine("%s/plots/kern_QED_num_reco.xmg",ens.path.c_str()),VV_QED,A_V,Z_V,SL_V,M_V,resc_a);
      cout<<"bare amu_QED: "<<bare_QED_correl.ave_err()<<" + "<<bare_QED_remaind.ave_err()<<" = "<<bare_QED[iens].ave_err()<<endl;
      
      ratio[iens]=bare_QED[iens]/bare_LO[iens];
      cout<<" Ratio: "<<ratio[iens].ave_err()<<endl;
    }
  
  vector<string> beta_list={"1.90","1.95","2.10"};
  
  //loop over analysis flags and input scale determination
  dbvec_t cLO(ninput_an*nan_syst);
  dbvec_t cQED(ninput_an*nan_syst);
  dbvec_t cRAT(ninput_an*nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	prepare_az(input_an_id);
	dboot_t &f0=lat_par[input_an_id].f0;
	dboot_t &B0=lat_par[input_an_id].B0;
	
	vector<cont_chir_fit_data_t> data_LO;
	vector<cont_chir_fit_data_t> data_QED;
	vector<cont_chir_fit_data_t> data_ratio;
	for(size_t iens=0;iens<ens_data.size();iens++)
	  {
	    ens_data_t &ens=ens_data[iens];
	    int ib=ens.ib;
	    
	    //set bootstrap
	    bi=jack_index[input_an_id][ens.iult];
	    
	    dboot_t LO=dboot_t(bi,bare_LO[iens])*sqr(Za[ib]);
	    data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,ens_data[iens].ib,ens_data[iens].L,LO,LO));
	    
	    dboot_t QED=dboot_t(bi,bare_QED[iens])*sqr(Za[ib]);
	    data_QED.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,ens_data[iens].ib,ens_data[iens].L,QED,QED));
	    
	    dboot_t RAT=dboot_t(bi,ratio[iens]);
	    data_ratio.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,ens_data[iens].ib,ens_data[iens].L,RAT,RAT));
	  }
	
	int iai=ind_an({input_an_id,an_flag});
	
	cLO[iai]=cont_chir_fit_LO(alist,zlist,f0,B0,data_LO,lat_par[input_an_id].ml,combine("plots/cont_chir_LO_flag%zu_an%zu.xmg",an_flag,input_an_id,"%s"),an_flag,false,beta_list);
	cout<<cLO[iai].ave_err()<<endl;
	//cQED[iai]=cont_chir_fit_QED(alist,zlist,data_QED,lat_par[input_an_id].ml,combine("plots/cont_chir_QED_flag%zu_an%zu.xmg",an_flag,input_an_id,"%s"),an_flag,false,beta_list);
	//cRAT[iai]=cont_chir_fit_ratio(alist,zlist,data_ratio,lat_par[input_an_id].ml,combine("plots/cont_chir_rat_flag%zu_an%zu.xmg",an_flag,input_an_id,"%s"),an_flag,false,beta_list);
      }

   perform_analysis(cLO,"LO");
   //cout<<crat.ave_err()<<endl;
   //cout<<cQED.ave_err()<<endl;
   
   return 0;
}
