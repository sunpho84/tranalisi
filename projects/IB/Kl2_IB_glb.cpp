#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <Kl2_IB_fit.hpp>
#include <common.hpp>

index_t<2> ind_an;

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

template <class T1,class T2> T1 FVE_M2(const T1 &M,const T2 &L)
{
  const double FVE_k=2.837297;
  return -FVE_k*alpha_em/L*(M+2.0/L);
}

template <class T1,class T2> T1 FVE_M(const T1 &M,const T2 &L)
{
  const double FVE_k=2.837297;
  return -FVE_k*alpha_em/L/2.0*(1.0+2.0/L/M);
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

double syst_analysis(const vector<ave_err_t> &v)
{
  ave_err_t ae;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave;
      ae.ave+=a;
      ae.err+=sqr(a);
    }
  ae.ave/=v.size();
  ae.err/=v.size();
  ae.err-=sqr(ae.ave);
  ae.err=sqrt(fabs(ae.err));
  
  return ae.err;
}

ave_err_t stat_analysis(const vector<ave_err_t> &v)
{
  ave_err_t ae;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave;
      double e=v[i].err;
      ae.ave+=a;
      ae.err+=sqr(e);
    }
  ae.ave/=v.size();
  ae.err/=v.size();
  ae.err=sqrt(fabs(ae.err));
  
  return ae;
}

vector<ave_err_t> ave_analyses(const dbvec_t &v)
{
  vector<ave_err_t> input_an_ave_err(nan_syst);
  dbvec_t v_an(ninput_an);
  
  for(size_t i=0;i<nan_syst;i++)
    {
      for(size_t j=0;j<ninput_an;j++)
	v_an[j]=v[j+i*ninput_an];
      input_an_ave_err[i]=eq_28_analysis(v_an);
    }
  
  return input_an_ave_err;
}

void syst_analysis_sep(const vector<ave_err_t> &v)
{
  double db[12];
  
  db[0]=(v[0].ave-v[1].ave)/2.0;
  db[1]=(v[2].ave-v[3].ave)/2.0;
  db[2]=(v[4].ave-v[5].ave)/2.0;
  db[3]=(v[6].ave-v[7].ave)/2.0;
  db[4]=(v[0].ave-v[2].ave)/2.0;
  db[5]=(v[1].ave-v[3].ave)/2.0;
  db[6]=(v[4].ave-v[6].ave)/2.0;
  db[7]=(v[5].ave-v[7].ave)/2.0;
  db[8]=(v[0].ave-v[4].ave)/2.0;
  db[9]=(v[1].ave-v[5].ave)/2.0;
  db[10]=(v[2].ave-v[6].ave)/2.0;
  db[11]=(v[3].ave-v[7].ave)/2.0;

  double S2cont,S2chir,S2fse;

  S2cont=1.0/24.0*(pow(db[0],2.0)+pow(db[1],2.0)+pow(db[2],2.0)+pow(db[3],2.0)+pow(db[0]+db[1],2.0)+pow(db[0]+db[2],2.0)+pow(db[1]+db[3],2.0)+pow(db[2]+db[3],2.0))+1.0/48.0*(pow(db[0]+db[3],2.0)+pow(db[1]+db[2],2.0));
  
  S2chir=1.0/24.0*(pow(db[8],2.0)+pow(db[9],2.0)+pow(db[10],2.0)+pow(db[11],2.0)+pow(db[8]+db[9],2.0)+pow(db[8]+db[10],2.0)+pow(db[9]+db[11],2.0)+pow(db[10]+db[11],2.0))+1.0/48.0*(pow(db[8]+db[11],2.0)+pow(db[9]+db[10],2.0));
  
  S2fse=1.0/24.0*(pow(db[4],2.0)+pow(db[5],2.0)+pow(db[6],2.0)+pow(db[7],2.0)+pow(db[4]+db[5],2.0)+pow(db[4]+db[6],2.0)+pow(db[5]+db[7],2.0)+pow(db[6]+db[7],2.0))+1.0/48.0*(pow(db[4]+db[7],2.0)+pow(db[5]+db[6],2.0));

  cout<<"cont: "<<sqrt(S2cont)<<endl;
  cout<<"chir: "<<sqrt(S2chir)<<endl;
  cout<<"fse: "<<sqrt(S2fse)<<endl;
  
}

int main(int narg,char **arg)
{
  ind_an.set_ranges({ninput_an,nan_syst});
  
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
  dbvec_t aDeltam_cr_u(raw_data.size());
  dbvec_t aDeltam_cr_d(raw_data.size());
  dbvec_t aMD(raw_data.size());
  dbvec_t aMD_sl_exch(raw_data.size());
  dbvec_t aMD_sl_selftad(raw_data.size());
  dbvec_t MD_sl_s(raw_data.size());
  dbvec_t MD_sl_p(raw_data.size());
  dbvec_t aMD_sl_selftad_revins(raw_data.size());
  dbvec_t MD_sl_p_revins(raw_data.size());

  //output
  dbvec_t output_dM2Pi(ninput_an*nan_syst);
  dbvec_t output_dM2K_QED(ninput_an*nan_syst);
  dbvec_t output_dM2K_QCD_over_minus_two_Deltamud(ninput_an*nan_syst);
  dbvec_t output_dM2D_QED(ninput_an*nan_syst);
  dbvec_t output_dM2D_QCD(ninput_an*nan_syst);
  dbvec_t output_dM2D_QCD_ind(ninput_an*nan_syst);
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
  vector<ave_err_t> v_ave_an_dM2D_QCD_ind(nan_syst);
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
	    
	    dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*lat_par[input_an_id].Z[ibeta]*(6.0*log(mu_MS*a)-22.596)/(32.0*sqr(M_PI)));
	    aDeltam_cr_u[iens]=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(eu);
	    aDeltam_cr_d[iens]=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(ed);
	    
	    dboot_t aMPi=dboot_t(bi,raw_data[iens].pi_mass);
	    dboot_t aMPi_sl_exch=dboot_t(bi,raw_data[iens].pi_SL_exch);
	    dboot_t aMPi_sl_selftad=dboot_t(bi,raw_data[iens].pi_SL_selftad);
	    dboot_t MPi_sl_p=dboot_t(bi,raw_data[iens].pi_SL_p);
	    a2M2Pi[iens]=aMPi*aMPi;
	    da2M2Pi[iens]=aMPi*sqr(eu-ed)*e2*aMPi_sl_exch;
	    FVE_da2M2Pi[iens]=FVE_M2(aMPi,raw_data[iens].L);
	    
	    dboot_t aMK=dboot_t(bi,raw_data[iens].k_mass);
	    dboot_t aMK_sl_exch=dboot_t(bi,raw_data[iens].k_SL_exch);
	    dboot_t aMK_sl_selftad=dboot_t(bi,raw_data[iens].k_SL_selftad);
	    dboot_t MK_sl_s=dboot_t(bi,raw_data[iens].k_SL_s);
	    dboot_t MK_sl_p=dboot_t(bi,raw_data[iens].k_SL_p);
	    dboot_t aMK_sl_selftad_revins=dboot_t(bi,raw_data[iens].k_SL_selftad_revins);
	    dboot_t MK_sl_p_revins=dboot_t(bi,raw_data[iens].k_SL_p_revins);
	    
	    dboot_t daMK_QED=
	      -2.0*ml[iens]*a*MK_sl_s/Z_QED
	      -(aDeltam_cr_u[iens]-aDeltam_cr_d[iens])*MK_sl_p
	      +(sqr(eu)-sqr(ed))*e2*(aMK_sl_exch-aMK_sl_selftad);
	    da2M2K_QED[iens]=daMK_QED*2.0*aMK;
	    FVE_da2M2K[iens]=FVE_M2(aMK,raw_data[iens].L);
	    daM2K_QCD_over_minus_two_Deltamud[iens]=2.0*aMK*MK_sl_s;
	    
	    aMD[iens]=dboot_t(bi,raw_data[iens].D_mass);
	    aMD_sl_exch[iens]=dboot_t(bi,raw_data[iens].D_SL_exch);
	    aMD_sl_selftad[iens]=dboot_t(bi,raw_data[iens].D_SL_selftad);
	    MD_sl_s[iens]=dboot_t(bi,raw_data[iens].D_SL_s);
	    MD_sl_p[iens]=dboot_t(bi,raw_data[iens].D_SL_p);
	    aMD_sl_selftad_revins[iens]=dboot_t(bi,raw_data[iens].D_SL_selftad_revins);
	    MD_sl_p_revins[iens]=dboot_t(bi,raw_data[iens].D_SL_p_revins);
	    MD[iens]=aMD[iens]/a;
	    
	    dboot_t daMD_QED=
	      2.0*ml[iens]*a*MD_sl_s[iens]/Z_QED
	      -(aDeltam_cr_d[iens]-aDeltam_cr_u[iens])*MD_sl_p[iens]
	      +(sqr(eu)-sqr(ed))*e2*aMD_sl_selftad[iens]+(eu-ed)*eu*e2*aMD_sl_exch[iens];
	    da2M2D_QED[iens]=daMD_QED*2.0*aMD[iens];
	    FVE_da2M2D[iens]=FVE_M2(aMD[iens],raw_data[iens].L);
	    
	    dboot_t aMDs=dboot_t(bi,raw_data[iens].Ds_mass);
	    MDs[iens]=aMDs/a;
	    
	    epsilon_gamma[iens]=(da2M2K_QED[iens]/da2M2Pi[iens])-1.0;
	    epsilon_gamma_minusFVE[iens]=((da2M2K_QED[iens]-FVE_da2M2K[iens])/(da2M2Pi[iens]-FVE_da2M2Pi[iens]))-1.0;
	    
	    dboot_t num_epsilon_Pi0=2.0*aMPi*(-(sqr(eu)+sqr(ed))*e2*(aMPi_sl_exch/2.0+aMPi_sl_selftad)-(aDeltam_cr_u[iens]+aDeltam_cr_d[iens])*MPi_sl_p);
	    epsilon_Pi0[iens]=num_epsilon_Pi0/da2M2Pi[iens];
	    epsilon_Pi0_minusFVE[iens]=num_epsilon_Pi0/(da2M2Pi[iens]-FVE_da2M2Pi[iens]);
	    
	    dboot_t num_epsilon_K0=2.0*aMK*(-sqr(ed)*e2*(aMK_sl_exch+aMK_sl_selftad_revins+aMK_sl_selftad)-aDeltam_cr_d[iens]*(MK_sl_p+MK_sl_p_revins));
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
	    data_dM2Pi.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						      da2M2Pi[iens]-FVE_da2M2Pi[iens],da2M2Pi[iens]));
	
	output_dM2Pi[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2Pi(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2Pi,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2Pi_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
		
	cout<<"                      QED K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_dM2K_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
							 da2M2K_QED[iens]-FVE_da2M2K[iens],da2M2K_QED[iens]));
	
	output_dM2K_QED[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2K_QED(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QED,lat_par[input_an_id].ml,lat_par[input_an_id].ms,combine("plots/cont_chir_fit_dM2K_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                      QCD K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QCD_over_minus_two_Deltamud;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2K_QCD_over_minus_two_Deltamud.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
									       daM2K_QCD_over_minus_two_Deltamud[iens],daM2K_QCD_over_minus_two_Deltamud[iens]));
	
	output_dM2K_QCD_over_minus_two_Deltamud[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2K_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QCD_over_minus_two_Deltamud,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2K_QCD_over_minus_two_Deltamud_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         M2Pi0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2Pi0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_M2Pi0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
							 a2M2Pi0g[iens],a2M2Pi0g[iens]));
	
	output_M2Pi0g[input_an_id+an_flag*ninput_an]=cont_chir_fit_M2Pi0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2Pi0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2Pi0g_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         M2K0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2K0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_M2K0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
							 a2M2K0g[iens],a2M2K0g[iens]));
	
	output_M2K0g[input_an_id+an_flag*ninput_an]=cont_chir_fit_M2K0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2K0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2K0g_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                     Epsilon "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_epsilon.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						      epsilon_gamma_minusFVE[iens],epsilon_gamma[iens]));
	
	output_epsilon[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_gamma_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                        EpsilonPi0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_Pi0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_epsilon_Pi0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
							  epsilon_Pi0_minusFVE[iens],epsilon_Pi0[iens]));
	
	output_epsilon_Pi0[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon_Pi0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_Pi0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_Pi0_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         EpsilonK0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_K0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_epsilon_K0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
							 epsilon_K0_minusFVE[iens],epsilon_K0[iens]));
	
	output_epsilon_K0[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon_K0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_K0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_K0_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         QED D "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2D_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2D_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						      da2M2D_QED[iens]-FVE_da2M2D[iens],da2M2D_QED[iens]));
	
	output_dM2D_QED[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_dM2D_QED,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),"$$[M^2_{D^+}-M^2_{D^0})_{QED} [GeV^2]",2.0,0.0,an_flag,1,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	plot_ens_data(combine("plots/MD_an%zu.xmg",input_an_id),ml,MD);
	plot_ens_data(combine("plots/MD_sl_s_an%zu.xmg",input_an_id),ml,MD_sl_s);
	// plot_ens_data(combine("plots/MDs_an%zu.xmg",input_an_id),ml,MDs);
		
	}
  
  ////////////////////output//////////////////
  v_ave_an_dM2Pi=ave_analyses(output_dM2Pi);
  cout<<"dM2Pi: "<<stat_analysis(v_ave_an_dM2Pi)<<" "<<syst_analysis(v_ave_an_dM2Pi)<<endl;
  syst_analysis_sep(v_ave_an_dM2Pi);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2Pi: "<<v_ave_an_dM2Pi[i]<<endl;
  
  v_ave_an_dM2K_QED=ave_analyses(output_dM2K_QED);
  cout<<"dM2K_QED: "<<stat_analysis(v_ave_an_dM2K_QED)<<" "<<syst_analysis(v_ave_an_dM2K_QED)<<endl;
  syst_analysis_sep(v_ave_an_dM2K_QED);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QED: "<<v_ave_an_dM2K_QED[i]<<endl;
  
  v_ave_an_dM2K_QCD_over_minus_two_Deltamud=ave_analyses(output_dM2K_QCD_over_minus_two_Deltamud);
  cout<<"dM2K_QCD_over_minus_two_Deltamud: "<<stat_analysis(v_ave_an_dM2K_QCD_over_minus_two_Deltamud)<<" "<<syst_analysis(v_ave_an_dM2K_QCD_over_minus_two_Deltamud)<<endl;
  syst_analysis_sep(v_ave_an_dM2K_QCD_over_minus_two_Deltamud);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QCD_over_minus_two_Deltamud: "<<v_ave_an_dM2K_QCD_over_minus_two_Deltamud[i]<<endl;
  
  v_ave_an_M2Pi0g=ave_analyses(output_M2Pi0g);
  cout<<"M2Pi0g: "<<stat_analysis(v_ave_an_M2Pi0g)<<" "<<syst_analysis(v_ave_an_M2Pi0g)<<endl;
  syst_analysis_sep(v_ave_an_M2Pi0g);
  for(size_t i=0;i<8;i++)
    cout<<"an_M2Pi0g: "<<v_ave_an_M2Pi0g[i]<<endl;
  
  v_ave_an_M2K0g=ave_analyses(output_M2K0g);
  cout<<"M2K0g: "<<stat_analysis(v_ave_an_M2K0g)<<" "<<syst_analysis(v_ave_an_M2K0g)<<endl;
  syst_analysis_sep(v_ave_an_M2K0g);
  for(size_t i=0;i<8;i++)
    cout<<"an_M2K0g: "<<v_ave_an_M2K0g[i]<<endl;
  
  v_ave_an_epsilon=ave_analyses(output_epsilon);
  cout<<"Epsilon: "<<stat_analysis(v_ave_an_epsilon)<<" "<<syst_analysis(v_ave_an_epsilon)<<endl;
  syst_analysis_sep(v_ave_an_epsilon);
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon: "<<v_ave_an_epsilon[i]<<endl;
  
  v_ave_an_epsilon_Pi0=ave_analyses(output_epsilon_Pi0);
  cout<<"Epsilon_Pi0: "<<stat_analysis(v_ave_an_epsilon_Pi0)<<" "<<syst_analysis(v_ave_an_epsilon_Pi0)<<endl;
  syst_analysis_sep(v_ave_an_epsilon_Pi0);
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon_Pi0: "<<v_ave_an_epsilon_Pi0[i]<<endl;
  
  v_ave_an_epsilon_K0=ave_analyses(output_epsilon_K0);
  cout<<"Epsilon_K0: "<<stat_analysis(v_ave_an_epsilon_K0)<<" "<<syst_analysis(v_ave_an_epsilon_K0)<<endl;
  syst_analysis_sep(v_ave_an_epsilon_K0);
  for(size_t i=0;i<8;i++)
    cout<<"an_epsilon_K0: "<<v_ave_an_epsilon_K0[i]<<endl;
  
  v_ave_an_dM2D_QED=ave_analyses(output_dM2D_QED);
  cout<<"dM2D_QED: "<<stat_analysis(v_ave_an_dM2D_QED)<<" "<<syst_analysis(v_ave_an_dM2D_QED)<<endl;
  syst_analysis_sep(v_ave_an_dM2D_QED);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2D_QED: "<<v_ave_an_dM2D_QED[i]<<endl;
  
  /////////////////epsilon///////////////
  dbvec_t epsilon_Pi0_ind(ninput_an*nan_syst);
  dbvec_t epsilon_K0_ind(ninput_an*nan_syst);
  dbvec_t epsilon_gamma_ind(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_Pi0_ind(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_K0_ind(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_gamma_ind(nan_syst);
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	epsilon_Pi0_ind[input_an_id+an_flag*ninput_an]=dboot_t(output_M2Pi0g[input_an_id+an_flag*ninput_an]/output_dM2Pi[input_an_id+an_flag*ninput_an]);
	epsilon_K0_ind[input_an_id+an_flag*ninput_an]=dboot_t(output_M2K0g[input_an_id+an_flag*ninput_an]/output_dM2Pi[input_an_id+an_flag*ninput_an]);
	epsilon_gamma_ind[input_an_id+an_flag*ninput_an]=dboot_t(output_dM2K_QED[input_an_id+an_flag*ninput_an]/output_dM2Pi[input_an_id+an_flag*ninput_an]-1.0);
      }

  v_ave_an_epsilon_Pi0_ind=ave_analyses(epsilon_Pi0_ind);
  cout<<"Epsilon_Pi0_ind: "<<stat_analysis(v_ave_an_epsilon_Pi0_ind)<<" "<<syst_analysis(v_ave_an_epsilon_Pi0_ind)<<endl;
  syst_analysis_sep(v_ave_an_epsilon_Pi0_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_Pi0_ind: "<<v_ave_an_epsilon_Pi0_ind[i]<<endl;
  
  v_ave_an_epsilon_K0_ind=ave_analyses(epsilon_K0_ind);
  cout<<"Epsilon_K0_ind: "<<stat_analysis(v_ave_an_epsilon_K0_ind)<<" "<<syst_analysis(v_ave_an_epsilon_K0_ind)<<endl;
  syst_analysis_sep(v_ave_an_epsilon_K0_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_K0_ind: "<<v_ave_an_epsilon_K0_ind[i]<<endl;
  
  v_ave_an_epsilon_gamma_ind=ave_analyses(epsilon_gamma_ind);
  cout<<"Epsilon_gamma_ind: "<<stat_analysis(v_ave_an_epsilon_gamma_ind)<<" "<<syst_analysis(v_ave_an_epsilon_gamma_ind)<<endl;
  syst_analysis_sep(v_ave_an_epsilon_gamma_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Epsilon_gamma_ind: "<<v_ave_an_epsilon_gamma_ind[i]<<endl;

  /////////////////pion QCD//////////////
  const double mpi0=134.9766;
  const double mpip=139.57018;

  double M_pion_QCD=mpi0-v_ave_an_epsilon_Pi0[7].ave*(mpip-mpi0);
  double M_pion_QCD_ind=mpi0-v_ave_an_epsilon_Pi0_ind[7].ave*(mpip-mpi0);

  cout<<"pion QCD: "<<M_pion_QCD<<endl;
  cout<<"pion QCD_ind: "<<M_pion_QCD_ind<<endl;

  ///////////////////Deltamud/////////////////
  dbvec_t dM2K_QED_ind(ninput_an*nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QED_ind[input_an_id+an_flag*ninput_an]=dboot_t((output_epsilon[input_an_id+an_flag*ninput_an]+1.0)*output_dM2Pi[input_an_id+an_flag*ninput_an]);
      }
  
  dbvec_t dM2K_QCD(ninput_an*nan_syst);
  dbvec_t dM2K_QCD_ind(ninput_an*nan_syst);
  dbvec_t Deltamud(ninput_an*nan_syst);
  dbvec_t Deltamud_ind(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud(nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud_ind(nan_syst);
  dboot_t dM2K_exp;
  for(size_t iboot=0;iboot<nboots;iboot++)
    dM2K_exp[iboot]=-3.903;
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QCD[input_an_id+an_flag*ninput_an]=dM2K_exp-dboot_t(output_dM2K_QED[input_an_id+an_flag*ninput_an]*1000.0);
	dM2K_QCD_ind[input_an_id+an_flag*ninput_an]=dM2K_exp-dboot_t(dM2K_QED_ind[input_an_id+an_flag*ninput_an]*1000.0);
	Deltamud[input_an_id+an_flag*ninput_an]=dboot_t(-dM2K_QCD[input_an_id+an_flag*ninput_an]/output_dM2K_QCD_over_minus_two_Deltamud[input_an_id+an_flag*ninput_an]/2.0);
	Deltamud_ind[input_an_id+an_flag*ninput_an]=dboot_t(-dM2K_QCD_ind[input_an_id+an_flag*ninput_an]/output_dM2K_QCD_over_minus_two_Deltamud[input_an_id+an_flag*ninput_an]/2.0);
      }
  
  v_ave_an_Deltamud=ave_analyses(Deltamud);
  cout<<"Deltamud: "<<stat_analysis(v_ave_an_Deltamud)<<" "<<syst_analysis(v_ave_an_Deltamud)<<endl;
  syst_analysis_sep(v_ave_an_Deltamud);
  for(size_t i=0;i<8;i++)
    cout<<"an_Deltamud: "<<v_ave_an_Deltamud[i]<<endl;

  v_ave_an_Deltamud_ind=ave_analyses(Deltamud_ind);
  cout<<"Deltamud_ind: "<<stat_analysis(v_ave_an_Deltamud_ind)<<" "<<syst_analysis(v_ave_an_Deltamud_ind)<<endl;
  syst_analysis_sep(v_ave_an_Deltamud_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Deltamud_ind: "<<v_ave_an_Deltamud_ind[i]<<endl;
  
  /////////////////R,Q2,mu/md/////////////////////
  dbvec_t R(ninput_an*nan_syst);
  dbvec_t R_ind(ninput_an*nan_syst);
  dbvec_t Q2(ninput_an*nan_syst);
  dbvec_t Q2_ind(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud_ind(ninput_an*nan_syst);
  dbvec_t ratio_mu_md(ninput_an*nan_syst);
  dbvec_t ratio_mu_md_ind(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_R(nan_syst);
  vector<ave_err_t> v_ave_an_R_ind(nan_syst);
  vector<ave_err_t> v_ave_an_Q2(nan_syst);
  vector<ave_err_t> v_ave_an_Q2_ind(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md_ind(nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	R[input_an_id+an_flag*ninput_an]=dboot_t((lat_par[input_an_id].ms-lat_par[input_an_id].ml)/(2.0*Deltamud[input_an_id+an_flag*ninput_an]*1.0e-3));
	R_ind[input_an_id+an_flag*ninput_an]=dboot_t((lat_par[input_an_id].ms-lat_par[input_an_id].ml)/(2.0*Deltamud_ind[input_an_id+an_flag*ninput_an]*1.0e-3));
	Q2[input_an_id+an_flag*ninput_an]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Deltamud[input_an_id+an_flag*ninput_an]*lat_par[input_an_id].ml*1.0e-3));
	Q2_ind[input_an_id+an_flag*ninput_an]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Deltamud_ind[input_an_id+an_flag*ninput_an]*lat_par[input_an_id].ml*1.0e-3));
	Deltamud_over_mud[input_an_id+an_flag*ninput_an]=dboot_t(Deltamud[input_an_id+an_flag*ninput_an]*1.0e-3/lat_par[input_an_id].ml);
	Deltamud_over_mud_ind[input_an_id+an_flag*ninput_an]=dboot_t(Deltamud_ind[input_an_id+an_flag*ninput_an]*1.0e-3/lat_par[input_an_id].ml);
	ratio_mu_md[input_an_id+an_flag*ninput_an]=dboot_t((1.0-Deltamud_over_mud[input_an_id+an_flag*ninput_an])/(1.0+Deltamud_over_mud[input_an_id+an_flag*ninput_an]));
	ratio_mu_md_ind[input_an_id+an_flag*ninput_an]=dboot_t((1.0-Deltamud_over_mud_ind[input_an_id+an_flag*ninput_an])/(1.0+Deltamud_over_mud_ind[input_an_id+an_flag*ninput_an]));
      }
  
  v_ave_an_R=ave_analyses(R);
  cout<<"R: "<<stat_analysis(v_ave_an_R)<<" "<<syst_analysis(v_ave_an_R)<<endl;
  syst_analysis_sep(v_ave_an_R);
  for(size_t i=0;i<8;i++)
    cout<<"an_R: "<<v_ave_an_R[i]<<endl;

  v_ave_an_Q2=ave_analyses(Q2);
  cout<<"Q2: "<<stat_analysis(v_ave_an_Q2)<<" "<<syst_analysis(v_ave_an_Q2)<<endl;
  syst_analysis_sep(v_ave_an_Q2);
  for(size_t i=0;i<8;i++)
    cout<<"an_Q2: "<<v_ave_an_Q2[i]<<endl;
  
  v_ave_an_ratio_mu_md=ave_analyses(ratio_mu_md);
  cout<<"mu/md: "<<stat_analysis(v_ave_an_ratio_mu_md)<<" "<<syst_analysis(v_ave_an_ratio_mu_md)<<endl;
  syst_analysis_sep(v_ave_an_ratio_mu_md);
  for(size_t i=0;i<8;i++)
    cout<<"an_ratio_mu_md: "<<v_ave_an_ratio_mu_md[i]<<endl;

  v_ave_an_R_ind=ave_analyses(R_ind);
  cout<<"R_ind: "<<stat_analysis(v_ave_an_R_ind)<<" "<<syst_analysis(v_ave_an_R_ind)<<endl;
  syst_analysis_sep(v_ave_an_R_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_R_ind: "<<v_ave_an_R_ind[i]<<endl;

  v_ave_an_Q2_ind=ave_analyses(Q2_ind);
  cout<<"Q2_ind: "<<stat_analysis(v_ave_an_Q2_ind)<<" "<<syst_analysis(v_ave_an_Q2_ind)<<endl;
  syst_analysis_sep(v_ave_an_Q2_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Q2_ind: "<<v_ave_an_Q2_ind[i]<<endl;

  v_ave_an_ratio_mu_md_ind=ave_analyses(ratio_mu_md_ind);
  cout<<"mu/md_ind: "<<stat_analysis(v_ave_an_ratio_mu_md_ind)<<" "<<syst_analysis(v_ave_an_ratio_mu_md_ind)<<endl;
  syst_analysis_sep(v_ave_an_ratio_mu_md_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_ratio_mu_md_ind: "<<v_ave_an_ratio_mu_md_ind[i]<<endl;

  ////////////////////D meson////////////////////
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t daM2D_QCD(raw_data.size());
	dbvec_t daM2D_QCD_ind(raw_data.size());
      
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    daM2D_QCD[iens]=4.0*Deltamud[input_an_id+an_flag*ninput_an]*MD_sl_s[iens]*aMD[iens]/1000.0;
	    daM2D_QCD_ind[iens]=4.0*Deltamud_ind[input_an_id+an_flag*ninput_an]*MD_sl_s[iens]*aMD[iens]/1000.0;
	  }

	cout<<"-----------------------------------------------"<<endl;
	cout<<"                        an_flag: "<<an_flag<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<"                    input_an_id: "<<input_an_id<<endl;
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         QCD D "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dM2D_QCD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2D_QCD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						     daM2D_QCD[iens],daM2D_QCD[iens]));
      
	output_dM2D_QCD[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_dM2D_QCD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QCD_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$[M^2_{D^+}-M^2_{D^0}]_{QCD} [GeV^2]",1.0,1.0,an_flag,0,cov_flag,beta_list);
      
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         QCD D_ind "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dM2D_QCD_ind;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2D_QCD_ind.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						     daM2D_QCD_ind[iens],daM2D_QCD_ind[iens]));
      
	output_dM2D_QCD_ind[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_dM2D_QCD_ind,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QCD_ind_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$[M^2_{D^+}-M^2_{D^0}]_{QCD} [GeV^2]",1.0,1.0,an_flag,0,cov_flag,beta_list);
      
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
      }

  v_ave_an_dM2D_QCD=ave_analyses(output_dM2D_QCD);
  cout<<"dM2D_QCD: "<<stat_analysis(v_ave_an_dM2D_QCD)<<" "<<syst_analysis(v_ave_an_dM2D_QCD)<<endl;
  syst_analysis_sep(v_ave_an_dM2D_QCD);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2D_QCD: "<<v_ave_an_dM2D_QCD[i]<<endl;

  v_ave_an_dM2D_QCD_ind=ave_analyses(output_dM2D_QCD_ind);
  cout<<"dM2D_QCD_ind: "<<stat_analysis(v_ave_an_dM2D_QCD_ind)<<" "<<syst_analysis(v_ave_an_dM2D_QCD_ind)<<endl;
  syst_analysis_sep(v_ave_an_dM2D_QCD_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2D_QCD_ind: "<<v_ave_an_dM2D_QCD_ind[i]<<endl;
  
  dbvec_t dM2D(ninput_an*nan_syst);
  dbvec_t dM2D_ind(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_dM2D(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_ind(nan_syst);

  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2D[input_an_id+an_flag*ninput_an]=output_dM2D_QCD[input_an_id+an_flag*ninput_an]+output_dM2D_QED[input_an_id+an_flag*ninput_an];
	dM2D_ind[input_an_id+an_flag*ninput_an]=output_dM2D_QCD_ind[input_an_id+an_flag*ninput_an]+output_dM2D_QED[input_an_id+an_flag*ninput_an];
      }
  
  v_ave_an_dM2D=ave_analyses(dM2D);
  cout<<"M^2_{D^+}-M^2_{D^0}: "<<stat_analysis(v_ave_an_dM2D)<<" "<<syst_analysis(v_ave_an_dM2D)<<endl;
  syst_analysis_sep(v_ave_an_dM2D);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2D: "<<v_ave_an_dM2D[i]<<endl;
  
  v_ave_an_dM2D_ind=ave_analyses(dM2D_ind);
  cout<<"M^2_{D^+}-M^2_{D^0}_ind: "<<stat_analysis(v_ave_an_dM2D_ind)<<" "<<syst_analysis(v_ave_an_dM2D_ind)<<endl;
  syst_analysis_sep(v_ave_an_dM2D_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2D_ind: "<<v_ave_an_dM2D_ind[i]<<endl;
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t saMD(raw_data.size());
	dbvec_t FVE_saMD(raw_data.size());
	
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	  saMD[iens]=-(eu+ed)*eu*e2*aMD_sl_exch[iens]-2.0*sqr(eu)*e2*aMD_sl_selftad_revins[iens]-(sqr(eu)+sqr(ed))*e2*aMD_sl_selftad[iens]-2.0*aDeltam_cr_u[iens]*MD_sl_p_revins[iens]-(aDeltam_cr_u[iens]+aDeltam_cr_d[iens])*MD_sl_p[iens];
	  FVE_saMD[iens]=FVE_M(aMD[iens],raw_data[iens].L);
	  }

	cout<<"-----------------------------------------------"<<endl;
	cout<<"                        an_flag: "<<an_flag<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	cout<<"-----------------------------------------------"<<endl;
	cout<<"                    input_an_id: "<<input_an_id<<endl;
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         SMD "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_sMD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_sMD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens],raw_data[iens].ibeta,raw_data[iens].L,
						     saMD[iens]-FVE_saMD[iens],saMD[iens]));
      
	output_sMD[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_sMD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_sMD_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$M_{D^+}+M_{D^0} [GeV]",1.0,0.0,an_flag,1,cov_flag,beta_list);
      
	cout<<"-----------------------------------------------"<<endl;
      }

  v_ave_an_sMD=ave_analyses(output_sMD);
  cout<<"M_{D^+}+M_{D^0}: "<<stat_analysis(v_ave_an_sMD)<<" "<<syst_analysis(v_ave_an_sMD)<<endl;
  syst_analysis_sep(v_ave_an_sMD);
  for(size_t i=0;i<8;i++)
    cout<<"an_sMD: "<<v_ave_an_sMD[i]<<endl;
  
  return 0;
}
