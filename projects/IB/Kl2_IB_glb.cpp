#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <Kl2_IB_fit.hpp>
#include <common.hpp>

class ens_data_t
{
public:
  size_t iult,ibeta,L,useforL;
  double aml,ams,amc;
  string path;
  
  djack_t pi_mass,pi_SL_exch,pi_SL_selftad,pi_SL_s,pi_SL_p,k_mass,k_SL_exch,k_SL_selftad,k_SL_selftad_revins,k_SL_s,k_SL_p,k_SL_p_revins,D_mass,D_SL_exch,D_SL_selftad,D_SL_selftad_revins,D_SL_s,D_SL_p,D_SL_p_revins,Ds_mass;
  djack_t deltam_cr;
};

vector<ens_data_t> raw_data;

template <class T1,class T2> T1 FVE_d2M(const T1 &M,const T2 &L)
{
  const double FVE_k=2.837297;
  return -FVE_k*alpha_em/L*(M+2.0/L);
}

template <class T1,class T2> T1 FVE_dM(const T1 &M,const T2 &L)
{
  const double FVE_k=2.837297;
  return -FVE_k*alpha_em/L/2.0*(1+2.0/L/M);
}

//! plot the three ensemble separately
template <class Tx,class Ty> void plot_ens_data(grace_file_t &file,const Tx &x,const Ty &y)
{
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      file.new_data_set();
      for(size_t iens=0;iens<raw_data.size();iens++)
	if(raw_data[iens].ibeta==ibeta)
	  file<<x[iens].ave()<<" "<<y[iens].ave_err()<<endl;
    }
}

template <class Tx,class Ty> void plot_ens_data(string path,const Tx &x,const Ty &y)
{
  grace_file_t file(path);
  plot_ens_data(file,x,y);
}

int main(int narg,char **arg)
{
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  //read where to read input and how many ensemble
  string ens_pars=input.read<string>("UltimatePath");
  init_common_IB(ens_pars);
  size_t nens_used=input.read<int>("NEnsamble");
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t temp;
      
      //read the parameters of the ensemble
      temp.iult=input.read<size_t>("Ens");
      temp.ibeta=input.read<size_t>("beta");
      temp.L=input.read<size_t>("L");
      temp.useforL=input.read<size_t>("useforL");
      temp.aml=input.read<double>("aml");
      temp.ams=input.read<double>("ams");
      temp.amc=input.read<double>("amc");
      temp.path=input.read<string>("path");
      
      //read the observable of the pion
      raw_file_t obs_pi_file(combine("%s/pi_obs",temp.path.c_str()),"r");
      obs_pi_file.bin_read(temp.pi_mass);
      obs_pi_file.bin_read(temp.pi_SL_exch);
      obs_pi_file.bin_read(temp.pi_SL_selftad);
      obs_pi_file.bin_read(temp.pi_SL_s);
      obs_pi_file.bin_read(temp.pi_SL_p);

      //read the observable of the kaon
      raw_file_t obs_k_file(combine("%s/k_obs",temp.path.c_str()),"r");
      obs_k_file.bin_read(temp.k_mass);
      obs_k_file.bin_read(temp.k_SL_exch);
      obs_k_file.bin_read(temp.k_SL_selftad);
      obs_k_file.bin_read(temp.k_SL_s);
      obs_k_file.bin_read(temp.k_SL_p);
      obs_k_file.bin_read(temp.k_SL_selftad_revins);
      obs_k_file.bin_read(temp.k_SL_p_revins);

      //read the observable of the D meson
      raw_file_t obs_D_file(combine("%s/D_obs",temp.path.c_str()),"r");
      obs_D_file.bin_read(temp.D_mass);
      obs_D_file.bin_read(temp.D_SL_exch);
      obs_D_file.bin_read(temp.D_SL_selftad);
      obs_D_file.bin_read(temp.D_SL_s);
      obs_D_file.bin_read(temp.D_SL_p);
      obs_D_file.bin_read(temp.D_SL_selftad_revins);
      obs_D_file.bin_read(temp.D_SL_p_revins);

      //read the observable of the Ds meson
      raw_file_t obs_Ds_file(combine("%s/Ds_obs",temp.path.c_str()),"r");
      obs_Ds_file.bin_read(temp.Ds_mass);

      //read deltam_cr (ud)
      raw_file_t deltam_cr_file(combine("%s/ud_fit_deltam_cr",temp.path.c_str()),"r");
      deltam_cr_file.bin_read(temp.deltam_cr);
      
      //store in the raw_data vector
      raw_data.push_back(temp);
    }

  //alist, zlist
  dbvec_t alist(nbeta),zlist(nbeta);
  const vector<string> beta_list={"1.90","1.95","2.10"};
  
  //D meson
  dboot_t aMD;
  dboot_t aMD_sl_exch,aMD_sl_selftad,MD_sl_s,MD_sl_p,aMD_sl_selftad_revins,MD_sl_p_revins;

  //output
  dbvec_t output_dM2Pi(ninput_an*nan_syst);
  dbvec_t output_dM2K_QED(ninput_an*nan_syst);
  dbvec_t output_dM2K_QCD_over_minus_two_Deltamud(ninput_an*nan_syst);
  dbvec_t output_dM2D_QED(ninput_an*nan_syst);
  dbvec_t output_dM2D_QCD(ninput_an*nan_syst);
  dbvec_t output_sMD(ninput_an*nan_syst);
  dbvec_t output_epsilon(ninput_an*nan_syst);
  dbvec_t output_epsilon_Pi0(ninput_an*nan_syst);
  dbvec_t output_epsilon_K0(ninput_an*nan_syst);
  dbvec_t output_M2Pi0g(ninput_an*nan_syst);
  dbvec_t output_M2K0g(ninput_an*nan_syst);

  vector<ave_err_t> v_ave_an_dM2Pi(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QCD_over_minus_two_Deltamud(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_QCD(nan_syst);
  vector<ave_err_t> v_ave_an_sMD(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_Pi0(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_K0(nan_syst);
  vector<ave_err_t> v_ave_an_M2Pi0g(nan_syst);
  vector<ave_err_t> v_ave_an_M2K0g(nan_syst);
  
  bool cov_flag=false;

  //loop over analysis flags and input scale determination
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t ml(raw_data.size());
	dbvec_t a2M2Pi(raw_data.size());
	dbvec_t da2M2Pi(raw_data.size());
	dbvec_t FVE_da2M2Pi(raw_data.size());
	dbvec_t da2M2K_QED(raw_data.size());
	dbvec_t FVE_da2M2K(raw_data.size());
	dbvec_t daM2K_QCD_over_minus_two_Deltamud(raw_data.size());
	dbvec_t da2M2D_QED(raw_data.size());
	dbvec_t FVE_da2M2D(raw_data.size());
	dbvec_t MD(raw_data.size());
	dbvec_t MDs(raw_data.size());
	dbvec_t epsilon_gamma(raw_data.size());
	dbvec_t epsilon_gamma_minusFVE(raw_data.size());
	dbvec_t epsilon_Pi0(raw_data.size());
	dbvec_t epsilon_Pi0_minusFVE(raw_data.size());
	dbvec_t epsilon_K0(raw_data.size());
	dbvec_t epsilon_K0_minusFVE(raw_data.size());
	dbvec_t a2M2Pi0g(raw_data.size());
	dbvec_t a2M2K0g(raw_data.size());
	
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    size_t ens_id=raw_data[iens].iult;
	    size_t ibeta=raw_data[iens].ibeta;
	    
	    dboot_t a=1.0/lat_par[input_an_id].ainv[ibeta];
	    dboot_t Lphys=raw_data[iens].L*a;
	    ml[iens]=raw_data[iens].aml/lat_par[input_an_id].Z[ibeta]/a;
	    
	    boot_init_t &bi=jack_index[input_an_id][ens_id];
	    
	    dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*lat_par[input_an_id].Z[ibeta]*(6*log(mu_MS*a)-22.596)/(32*sqr(M_PI)));
	    dboot_t aDeltam_cr_u=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(eu);
	    dboot_t aDeltam_cr_d=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(ed);
	    
	    dboot_t aMPi=dboot_t(bi,raw_data[iens].pi_mass);
	    dboot_t aMPi_sl_exch=dboot_t(bi,raw_data[iens].pi_SL_exch);
	    dboot_t aMPi_sl_selftad=dboot_t(bi,raw_data[iens].pi_SL_selftad);
	    dboot_t MPi_sl_p=dboot_t(bi,raw_data[iens].pi_SL_p);
	    a2M2Pi[iens]=aMPi*aMPi;
	    da2M2Pi[iens]=aMPi*sqr(eu-ed)*e2*aMPi_sl_exch;
	    FVE_da2M2Pi[iens]=FVE_d2M(aMPi,raw_data[iens].L);
	    
	    dboot_t aMK=dboot_t(bi,raw_data[iens].k_mass);
	    dboot_t aMK_sl_exch=dboot_t(bi,raw_data[iens].k_SL_exch);
	    dboot_t aMK_sl_selftad=dboot_t(bi,raw_data[iens].k_SL_selftad);
	    dboot_t MK_sl_s=dboot_t(bi,raw_data[iens].k_SL_s);
	    dboot_t MK_sl_p=dboot_t(bi,raw_data[iens].k_SL_p);
	    dboot_t aMK_sl_selftad_revins=dboot_t(bi,raw_data[iens].k_SL_selftad_revins);
	    dboot_t MK_sl_p_revins=dboot_t(bi,raw_data[iens].k_SL_p_revins);
	    
	    dboot_t daMK_QED=
	      -2*ml[iens]*a*MK_sl_s/Z_QED
	      -(aDeltam_cr_u-aDeltam_cr_d)*MK_sl_p
	      +(sqr(eu)-sqr(ed))*e2*(aMK_sl_exch-aMK_sl_selftad);
	    da2M2K_QED[iens]=daMK_QED*2*aMK;
	    FVE_da2M2K[iens]=FVE_d2M(aMK,raw_data[iens].L);
	    daM2K_QCD_over_minus_two_Deltamud[iens]=2*aMK*MK_sl_s;
	    
	    aMD=dboot_t(bi,raw_data[iens].D_mass);
	    aMD_sl_exch=dboot_t(bi,raw_data[iens].D_SL_exch);
	    aMD_sl_selftad=dboot_t(bi,raw_data[iens].D_SL_selftad);
	    MD_sl_s=dboot_t(bi,raw_data[iens].D_SL_s);
	    MD_sl_p=dboot_t(bi,raw_data[iens].D_SL_p);
	    aMD_sl_selftad_revins=dboot_t(bi,raw_data[iens].D_SL_selftad_revins);
	    MD_sl_p_revins=dboot_t(bi,raw_data[iens].D_SL_p_revins);
	    MD[iens]=aMD/a;
	    
	    dboot_t daMD_QED=
	      2*ml[iens]*a*MD_sl_s/Z_QED
	      -(aDeltam_cr_d-aDeltam_cr_u)*MD_sl_p
	      +(sqr(eu)-sqr(ed))*e2*aMD_sl_selftad+(eu-ed)*eu*e2*aMD_sl_exch;
	    da2M2D_QED[iens]=daMD_QED*2*aMD;
	    FVE_da2M2D[iens]=FVE_d2M(aMD,raw_data[iens].L);
	    
	    dboot_t aMDs=dboot_t(bi,raw_data[iens].Ds_mass);
	    MDs[iens]=aMDs/a;
	    
	    epsilon_gamma[iens]=(da2M2K_QED[iens]/da2M2Pi[iens])-1.0;
	    epsilon_gamma_minusFVE[iens]=((da2M2K_QED[iens]-FVE_da2M2K[iens])/(da2M2Pi[iens]-FVE_da2M2Pi[iens]))-1.0;
	    
	    dboot_t num_epsilon_Pi0=2*aMPi*(-(sqr(eu)+sqr(ed))*e2*(aMPi_sl_exch/2.0+aMPi_sl_selftad)-(aDeltam_cr_u+aDeltam_cr_d)*MPi_sl_p);
	    epsilon_Pi0[iens]=num_epsilon_Pi0/da2M2Pi[iens];
	    epsilon_Pi0_minusFVE[iens]=num_epsilon_Pi0/(da2M2Pi[iens]-FVE_da2M2Pi[iens]);
	    
	    dboot_t num_epsilon_K0=2*aMK*(-sqr(ed)*e2*(aMK_sl_exch+aMK_sl_selftad_revins+aMK_sl_selftad)-aDeltam_cr_d*(MK_sl_p+MK_sl_p_revins));
	    epsilon_K0[iens]=num_epsilon_K0/da2M2Pi[iens];
	    epsilon_K0_minusFVE[iens]=num_epsilon_K0/(da2M2Pi[iens]-FVE_da2M2Pi[iens]);
	    a2M2Pi0g[iens]=num_epsilon_Pi0;
	    a2M2K0g[iens]=num_epsilon_K0;
	  }
	
	//prepare the list of a and z
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  {
	    alist[ibeta]=1.0/lat_par[input_an_id].ainv[ibeta];
	    zlist[ibeta]=lat_par[input_an_id].Z[ibeta];
	  }

	cout<<"-----------------------------------------------"<<endl;
	cout<<"                        an_flag: "<<an_flag<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<"                    input_an_id: "<<input_an_id<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         Pi "<<endl;
	cout<<endl;
	
	//data to fit
	vector<cont_chir_fit_data_t> data_dM2Pi;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_dM2Pi.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						      da2M2Pi[iens]-FVE_da2M2Pi[iens],da2M2Pi[iens]));
	
	output_dM2Pi[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2Pi(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2Pi,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2Pi_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	/*	
	cout<<"                      QED K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_dM2K_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							 da2M2K_QED[iens]-FVE_da2M2K[iens],da2M2K_QED[iens]));
	
	output_dM2K_QED[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2K_QED(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QED,lat_par[input_an_id].ml,lat_par[input_an_id].ms,combine("plots/cont_chir_fit_dM2K_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                      QCD K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QCD_over_minus_two_Deltamud;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2K_QCD_over_minus_two_Deltamud.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
									       daM2K_QCD_over_minus_two_Deltamud[iens],daM2K_QCD_over_minus_two_Deltamud[iens]));
	
	output_dM2K_QCD_over_minus_two_Deltamud[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2K_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QCD_over_minus_two_Deltamud,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2K_QCD_over_minus_two_Deltamud_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         M2Pi0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2Pi0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_M2Pi0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							 a2M2Pi0g[iens],a2M2Pi0g[iens]));
	output_M2Pi0g[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2Pi0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2Pi0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2Pi0g_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         M2K0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2K0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_M2K0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							 a2M2K0g[iens],a2M2K0g[iens]));
	output_M2K0g[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2K0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2K0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2K0g_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
        cout<<"                     Epsilon "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						      epsilon_gamma_minusFVE[iens],epsilon_gamma[iens]));
	output_epsilon[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_gamma_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                        EpsilonPi0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_Pi0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon_Pi0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							  epsilon_Pi0_minusFVE[iens],epsilon_Pi0[iens]));
	output_epsilon_Pi0[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon_Pi0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_Pi0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_Pi0_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         EpsilonK0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_K0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon_K0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							 epsilon_K0_minusFVE[iens],epsilon_K0[iens]));
	output_epsilon_K0[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon_K0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_K0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_K0_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
		
	cout<<"                         QED D "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2D_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2D_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						      daMD_QED[iens],daMD_QED[iens]));
	
	output_dMD_QED[ind_an({input_an_id,an_flag})]=cont_chir_constant_fit(alist,zlist,data_dMD_QED,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dMD_QED_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$(M_{D^+}-M_{D^0})^{QED}[GeV]",1.0,0.0,an_flag,0,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	// plot_ens_data(combine("plots/MD_an%zu.xmg",input_an_id),ml,MD);
	// plot_ens_data(combine("plots/MDs_an%zu.xmg",input_an_id),ml,MDs);
	*/	
	}
  /*
  ////////////////////output//////////////////
  
  v_ave_an_dM2Pi=ave_analyses(output_dM2Pi);
  cout<<"dM2Pi: "<<stat_analysis(v_ave_an_dM2Pi)<<" "<<syst_analysis(v_ave_an_dM2Pi)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2Pi: "<<v_ave_an_dM2Pi[i]<<endl;
  
  v_ave_an_dM2K_QED=ave_analyses(output_dM2K_QED);
  cout<<"dM2K_QED: "<<stat_analysis(v_ave_an_dM2K_QED)<<" "<<syst_analysis(v_ave_an_dM2K_QED)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QED: "<<v_ave_an_dM2K_QED[i]<<endl;
  
  v_ave_an_dM2K_QCD_over_minus_two_Deltamud=ave_analyses(output_dM2K_QCD_over_minus_two_Deltamud);
  cout<<"dM2K_QCD_over_minus_two_Deltamud: "<<stat_analysis(v_ave_an_dM2K_QCD_over_minus_two_Deltamud)<<" "<<syst_analysis(v_ave_an_dM2K_QCD_over_minus_two_Deltamud)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QCD_over_minus_two_Deltamud: "<<v_ave_an_dM2K_QCD_over_minus_two_Deltamud[i]<<endl;
  
  v_ave_an_M2Pi0g=ave_analyses(output_M2Pi0g);
  cout<<"M2Pi0g: "<<stat_analysis(v_ave_an_M2Pi0g)<<" "<<syst_analysis(v_ave_an_M2Pi0g)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_M2Pi0g: "<<v_ave_an_M2Pi0g[i]<<endl;
  
  v_ave_an_M2K0g=ave_analyses(output_M2K0g);
  cout<<"M2K0g: "<<stat_analysis(v_ave_an_M2K0g)<<" "<<syst_analysis(v_ave_an_M2K0g)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_M2K0g: "<<v_ave_an_M2K0g[i]<<endl;
  
  v_ave_an_epsilon=ave_analyses(output_epsilon);
  cout<<"Epsilon: "<<stat_analysis(v_ave_an_epsilon)<<" "<<syst_analysis(v_ave_an_epsilon)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon: "<<v_ave_an_epsilon[i]<<endl;
  
  v_ave_an_epsilon_Pi0=ave_analyses(output_epsilon_Pi0);
  cout<<"Epsilon_Pi0: "<<stat_analysis(v_ave_an_epsilon_Pi0)<<" "<<syst_analysis(v_ave_an_epsilon_Pi0)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon_Pi0: "<<v_ave_an_epsilon_Pi0[i]<<endl;
  
  v_ave_an_epsilon_K0=ave_analyses(output_epsilon_K0);
  cout<<"Epsilon_K0: "<<stat_analysis(v_ave_an_epsilon_K0)<<" "<<syst_analysis(v_ave_an_epsilon_K0)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon_K0: "<<v_ave_an_epsilon_K0[i]<<endl;
  
  v_ave_an_dMD_QED=ave_analyses(output_dMD_QED);
  cout<<"dMD_QED: "<<stat_analysis(v_ave_an_dMD_QED)<<" "<<syst_analysis(v_ave_an_dMD_QED)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dMD_QED: "<<v_ave_an_dMD_QED[i]<<endl;
  
  /////////////////pion QCD//////////////
  const double mpi0=134.9766;
  const double mpip=139.57018;

  double M_pion_QCD=mpi0-v_ave_an_epsilon_Pi0[7].ave*(mpip-mpi0);

  cout<<"pion QCD: "<<M_pion_QCD<<endl;
  
  dbvec_t epsilon_Pi0_lat(ninput_an*nan_syst);
  dbvec_t epsilon_K0_lat(ninput_an*nan_syst);
  dbvec_t epsilon_gamma_lat(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_Pi0_lat(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_K0_lat(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_gamma_lat(nan_syst);
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	epsilon_Pi0_lat[ind_an({input_an_id,an_flag})]=dboot_t(output_M2Pi0g[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]);
	epsilon_K0_lat[ind_an({input_an_id,an_flag})]=dboot_t(output_M2K0g[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]);
	epsilon_gamma_lat[ind_an({input_an_id,an_flag})]=dboot_t(output_dM2K_QED[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]-1.0);
      }

  v_ave_an_epsilon_Pi0_lat=ave_analyses(epsilon_Pi0_lat);
  cout<<"Epsilon_Pi0_lat: "<<stat_analysis(v_ave_an_epsilon_Pi0_lat)<<" "<<syst_analysis(v_ave_an_epsilon_Pi0_lat)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_Pi0_lat: "<<v_ave_an_epsilon_Pi0_lat[i]<<endl;
  
  v_ave_an_epsilon_K0_lat=ave_analyses(epsilon_K0_lat);
  cout<<"Epsilon_K0_lat: "<<stat_analysis(v_ave_an_epsilon_K0_lat)<<" "<<syst_analysis(v_ave_an_epsilon_K0_lat)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_K0_lat: "<<v_ave_an_epsilon_K0_lat[i]<<endl;
  
  v_ave_an_epsilon_gamma_lat=ave_analyses(epsilon_gamma_lat);
  cout<<"Epsilon_gamma_lat: "<<stat_analysis(v_ave_an_epsilon_gamma_lat)<<" "<<syst_analysis(v_ave_an_epsilon_gamma_lat)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_gamma_lat: "<<v_ave_an_epsilon_gamma_lat[i]<<endl;
  
  dbvec_t dM2K_QED_lat(ninput_an*nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QED_lat[ind_an({input_an_id,an_flag})]=dboot_t((output_epsilon[ind_an({input_an_id,an_flag})]+1.0)*output_dM2Pi[ind_an({input_an_id,an_flag})]);
      }
  
  dbvec_t dM2K_QCD(ninput_an*nan_syst);
  dbvec_t Deltamud(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud(nan_syst);
  dboot_t dM2K_exp;
  for(size_t iboot=0;iboot<nboots;iboot++)
    dM2K_exp[iboot]=-3.903;
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QCD[ind_an({input_an_id,an_flag})]=dM2K_exp-dboot_t(output_dM2K_QED[ind_an({input_an_id,an_flag})]*1000);
	Deltamud[ind_an({input_an_id,an_flag})]=dboot_t(-dM2K_QCD[ind_an({input_an_id,an_flag})]/output_dM2K_QCD_over_minus_two_Deltamud[ind_an({input_an_id,an_flag})]/2.0);
      }
  
  /////////////////Deltamud///////////////////
  v_ave_an_Deltamud=ave_analyses(Deltamud);
  cout<<"Deltamud: "<<stat_analysis(v_ave_an_Deltamud)<<" "<<syst_analysis(v_ave_an_Deltamud)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_Deltamud: "<<v_ave_an_Deltamud[i]<<endl;
  
  /////////////////R,Q2,mu/md/////////////////////
  dbvec_t R(ninput_an*nan_syst);
  dbvec_t Q2(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud(ninput_an*nan_syst);
  dbvec_t ratio_mu_md(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_R(nan_syst);
  vector<ave_err_t> v_ave_an_Q2(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md(nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	R[ind_an({input_an_id,an_flag})]=dboot_t((lat_par[input_an_id].ms-lat_par[input_an_id].ml)/(2*Deltamud[ind_an({input_an_id,an_flag})]*1.0e-3));
	Q2[ind_an({input_an_id,an_flag})]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4*Deltamud[ind_an({input_an_id,an_flag})]*lat_par[input_an_id].ml*1.0e-3));
	Deltamud_over_mud[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud[ind_an({input_an_id,an_flag})]*1.0e-3/lat_par[input_an_id].ml);
	ratio_mu_md[ind_an({input_an_id,an_flag})]=dboot_t((1-Deltamud_over_mud[ind_an({input_an_id,an_flag})])/(1+Deltamud_over_mud[ind_an({input_an_id,an_flag})]));
      }
  v_ave_an_R=ave_analyses(R);
  cout<<"R: "<<stat_analysis(v_ave_an_R)<<" "<<syst_analysis(v_ave_an_R)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_R: "<<v_ave_an_R[i]<<endl;
  v_ave_an_Q2=ave_analyses(Q2);
  cout<<"Q2: "<<stat_analysis(v_ave_an_Q2)<<" "<<syst_analysis(v_ave_an_Q2)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_Q2: "<<v_ave_an_Q2[i]<<endl;
  v_ave_an_ratio_mu_md=ave_analyses(ratio_mu_md);
  cout<<"mu/md: "<<stat_analysis(v_ave_an_ratio_mu_md)<<" "<<syst_analysis(v_ave_an_ratio_mu_md)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_ratio_mu_md: "<<v_ave_an_ratio_mu_md[i]<<endl;
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t dMD_QCD(raw_data.size());
      
	for(size_t iens=0;iens<raw_data.size();iens++)
	  dMD_QCD[iens]=2*Deltamud[ind_an({input_an_id,an_flag})]*MD_sl_s/1000;

	cout<<"-----------------------------------------------"<<endl;
	cout<<"                        an_flag: "<<an_flag<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<"                    input_an_id: "<<input_an_id<<endl;
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         QCD D "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dMD_QCD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dMD_QCD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						     dMD_QCD[iens],dMD_QCD[iens]));
      
	output_dMD_QCD[ind_an({input_an_id,an_flag})]=cont_chir_constant_fit(alist,zlist,data_dMD_QCD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_constant_dMD_QCD_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$(M_{D^+}-M_{D^0})^{QCD}[GeV]",0.0,1.0,an_flag,0,cov_flag,beta_list);
      
	cout<<"-----------------------------------------------"<<endl;
      }

  //////////////////output/////////////////////
  v_ave_an_dMD_QCD=ave_analyses(output_dMD_QCD);
  cout<<"dMD_QCD: "<<stat_analysis(v_ave_an_dMD_QCD)<<" "<<syst_analysis(v_ave_an_dMD_QCD)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dMD_QCD: "<<v_ave_an_dMD_QCD[i]<<endl;
  
  //////////////////D meson//////////////////
  dbvec_t dMD(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_dMD(nan_syst);

  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      dMD[ind_an({input_an_id,an_flag})]=output_dMD_QCD[ind_an({input_an_id,an_flag})]+output_dMD_QED[ind_an({input_an_id,an_flag})];
  
  v_ave_an_dMD=ave_analyses(dMD);
  cout<<"M_{D^+}-M_{D^0}: "<<stat_analysis(v_ave_an_dMD)<<" "<<syst_analysis(v_ave_an_dMD)<<endl;
  for(size_t i=0;i<8;i++)
    cout<<"an_dMD: "<<v_ave_an_dMD[i]<<endl;
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t saMD(raw_data.size());
      
	for(size_t iens=0;iens<raw_data.size();iens++)
	  saMD[iens]=-(eu+ed)*eu*e2*aMD_sl_exch-2*sqr(eu)*e2*aMD_sl_selftad_revins-(sqr(eu)+sqr(ed))*e2*aMD_sl_selftad+2*aDeltam_cr_u*MD_sl_p_revins;

	cout<<"-----------------------------------------------"<<endl;
	cout<<"                        an_flag: "<<an_flag<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<"                    input_an_id: "<<input_an_id<<endl;
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         QCD D "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dMD_QCD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dMD_QCD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						     dMD_QCD[iens],dMD_QCD[iens]));
      
	output_dMD_QCD[ind_an({input_an_id,an_flag})]=cont_chir_constant_fit(alist,zlist,data_dMD_QCD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_constant_dMD_QCD_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$(M_{D^+}-M_{D^0})^{QCD}[GeV]",0.0,1.0,an_flag,0,cov_flag);
      
	cout<<"-----------------------------------------------"<<endl;
      }
	*/
  return 0;
}
