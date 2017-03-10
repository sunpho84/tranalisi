#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <gm2_IB_common.hpp>

const int DT=4; //!< number of points to eliminate from corr
const int im=istrange;

//! ansatz fit
template <class Tpars,class Tm,class Ta>
Tpars cont_chir_ansatz(const Tpars &C,const Tpars &Kpi,const Tm &ml,const Ta &a,const Tpars &adep,const Tpars &adep_ml,double L,const Tpars &L2dep,const size_t an_flag)
{return C+Kpi*ml+a*a*(adep+ml*adep_ml)+L2dep/sqr(L);}

//! perform the fit to the continuum limit
dboot_t cont_chir_fit(const dbvec_t &a,const dbvec_t &z,const vector<cont_chir_fit_data_t> &ext_data,const dboot_t &ml_phys,const string &path,size_t an_flag,bool cov_flag,const vector<string> &beta_list)
{
  //set_printlevel(3);
  
  boot_fit_t boot_fit;
  size_t nbeta=a.size();
  cont_chir_fit_pars_t pars(nbeta);
  
  ave_err_t adep_guess({6e-4,1e-4});
  //ave_err_t adep_guess({0,0.1});
  ave_err_t adep_ml_guess({0.0,0.0});
  //ave_err_t adep_ml_guess({0.0,0.1});
  pars.add_az_pars(a,z,boot_fit);
  
  pars.add_adep_pars(adep_guess,adep_ml_guess,boot_fit);
  pars.iC=boot_fit.add_fit_par(pars.C,"C",1e-9,1e-9);
  //pars.iC=boot_fit.add_fit_par(pars.C,"C",1,1e-3);
  pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",0,1e-9);
  //pars.iKPi=boot_fit.add_fit_par(pars.KPi,"KPi",0,0.1);
  pars.iL2dep=boot_fit.add_fit_par(pars.L2dep,"L2dep",0.0,0.07);
  
  //boot_fit.fix_par(pars.iadep);
  boot_fit.fix_par(pars.iL2dep);
  //boot_fit.fix_par(pars.iKPi);
  boot_fit.fix_par(pars.iadep_ml);
  
  cont_chir_fit_minimize(ext_data,pars,boot_fit,0.0,0.0,[an_flag](const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double ac,double L)
			 {return cont_chir_ansatz(p[pars.iC],p[pars.iKPi],ml,ac,p[pars.iadep],p[pars.iadep_ml],L,p[pars.iL2dep],an_flag);},cov_flag);
  
  double a_cont=0;
  dboot_t phys_res=cont_chir_ansatz(pars.C,pars.KPi,ml_phys,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L2dep,an_flag);
  cout<<"result: "<<phys_res.ave_err()<<endl;
  
  plot_chir_fit(path,ext_data,pars,
		[&pars,an_flag]
		(double x,size_t ib)
		{return cont_chir_ansatz<double,double,double>
		    (pars.C.ave(),pars.KPi.ave(),x,pars.fit_a[ib].ave(),pars.adep.ave(),pars.adep_ml.ave(),inf_vol,pars.L2dep.ave(),an_flag);},
		bind(cont_chir_ansatz<dboot_t,double,double>,pars.C,pars.KPi,_1,a_cont,pars.adep,pars.adep_ml,inf_vol,pars.L2dep,an_flag),
		[&ext_data,&pars]
		(size_t idata,bool without_with_fse,size_t ib)
		{return dboot_t(without_with_fse?ext_data[idata].wfse:ext_data[idata].wofse-without_with_fse*pars.L2dep/sqr(ext_data[idata].L));},
		ml_phys,phys_res,"$$a_\\mu",beta_list);
  
  return phys_res;
}

int main(int narg,char **arg)
{
  gm2_initialize(narg,arg);
  
  djvec_t LO(ens_data.size());
  djvec_t QED(ens_data.size());
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
      cout<<"Upto: "<<upto;
      djack_t LO_correl=integrate_corr_times_kern_up_to(VV_LO,ens.T,resc_a,im,upto);
      cout<<", after: "<<upto<<endl;
      djack_t LO_remaind=integrate_LO_reco_from(Z_V,M_V,resc_a,im,upto);
      LO[iens]=LO_correl+LO_remaind;
      compare_LO_num_reco(combine("%s/plots/kern_LO_num_reco.xmg",ens.path.c_str()),VV_LO,Z_V,M_V,resc_a);
      cout<<"amu: "<< LO_correl.ave_err()<<" + "<<LO_remaind.ave_err()<<" = "<<LO[iens].ave_err()<<endl;
      
      djack_t QED_correl=integrate_corr_times_kern_up_to(VV_QED,ens.T,resc_a,im,upto);
      djack_t QED_remaind=integrate_QED_reco_from(A_V,Z_V,SL_V,M_V,resc_a,im,upto);
      compare_QED_num_reco(combine("%s/plots/kern_QED_num_reco.xmg",ens.path.c_str()),VV_QED,A_V,Z_V,SL_V,M_V,resc_a);
      cout<<"amu_QED: "<<QED_correl.ave_err()<<" + "<<QED_remaind.ave_err()<<endl;
      
      // ratio[iens]=QED[iens]/LO[iens];
      // cout<<" Ratio: "<<ratio[iens].ave_err()<<endl;
    }
  
  // int input_an_id=0;
  // prepare_az(input_an_id);
  
  // const size_t an_flag=0;
  // vector<cont_chir_fit_data_t> data_LO;
  // for(size_t iens=0;iens<ens_data.size();iens++)
  //   {
  //     //set bootstrap
  //     bi=jack_index[input_an_id][ens.iult];
  //     data_LO.push_back(cont_chir_fit_data_t(ens_data[iens].aml,ens_data[iens].aml,ens_data[iens].ib,ens_data[iens].L,ratio[iens],ratio[iens]));
  //   }
  // cont_chir_fit(alist,zlist,data_LO,lat_par[input_an_id].ml,combine("plots/cont_chir_LO_flag%zu_an%zu.xmg",an_flag,input_an_id,"%s"),an_flag,false);
  
  return 0;
}
