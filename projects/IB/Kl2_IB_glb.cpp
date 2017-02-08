#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <Kl2_IB_fit.hpp>

size_t ninput_an=8,nbeta=3,nboots=100;

index_t<2> ind_an;

const double eu=2.0/3;
const double ed=-1.0/3;

const vector<size_t> symbol={grace::SQUARE,grace::CIRCLE,grace::DIAMOND};
const vector<size_t> color={grace::GREEN4,grace::RED,grace::BLUE};

class lat_par_t
{
public:
  dboot_t ml,ms,mc,r0,f0,B0;
  dbvec_t ainv,Z;
 
  lat_par_t() : ainv(nbeta),Z(nbeta) {}
};

dboot_t read_boot(const raw_file_t &file)
{
  dboot_t out;
  for(size_t ib=0;ib<nboots;ib++) file.read(out[ib]);
  return out;
}

class ens_data_t
{
public:
  size_t iult,ibeta,L;
  double aml,ams,amc;
  string path;
  
  djack_t pi_mass,pi_SL_exch,pi_SL_selftad,pi_SL_s,pi_SL_p,k_mass,k_SL_exch,k_SL_selftad,k_SL_selftad_revins,k_SL_s,k_SL_p,k_SL_p_revins,D_mass,D_SL_exch,D_SL_selftad,D_SL_s,D_SL_p,Ds_mass;
  djack_t deltam_cr;
};

vector<lat_par_t> lat_par(ninput_an);
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

//! perform the analysis according to eq.28
ave_err_t eq_28_analysis(const dbvec_t &v)
{
  ave_err_t ae;
  double sigma=0;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave();
      double e=v[i].err();
      ae.ave+=a;
      ae.err+=sqr(a);
      sigma+=sqr(e);
    }
  ae.ave/=v.size();
  ae.err/=v.size();
  sigma/=v.size();
  ae.err-=sqr(ae.ave);
  ae.err=sqrt(fabs(ae.err)+sigma);
  
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

int main(int narg,char **arg)
{
  set_njacks(15);
  ind_an.set_ranges({ninput_an,nan_syst});
  
  //open input file
  string name="input_global.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  //read where to read input and how many ensemble
  string ens_pars=input.read<string>("UltimatePath");
  size_t nens_used=input.read<int>("NEnsamble");
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      ens_data_t temp;
      
      //read the parameters of the ensemble
      temp.iult=input.read<size_t>("Ens");
      temp.ibeta=input.read<size_t>("beta");
      temp.L=input.read<size_t>("L");
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

      //read the observable of the Ds meson
      raw_file_t obs_Ds_file(combine("%s/Ds_obs",temp.path.c_str()),"r");
      obs_Ds_file.bin_read(temp.Ds_mass);

      //read deltam_cr (ud)
      raw_file_t deltam_cr_file(combine("%s/ud_fit_deltam_cr",temp.path.c_str()),"r");
      deltam_cr_file.bin_read(temp.deltam_cr);
      
      //store in the raw_data vector
      raw_data.push_back(temp);
    }
  
  //test the inputx
  //cout<<raw_data[0].pi_mass.ave_err()<<endl;

  raw_file_t file(ens_pars,"r");

  const int nens_total=15;
  boot_init_t jack_index[ninput_an][nens_total];

  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[input_an_id].ainv[ibeta][iboot]);
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an/2;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[input_an_id].Z[ibeta][iboot]);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=ninput_an/2;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[input_an_id].Z[ibeta][iboot]);
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_total;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[input_an_id][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].B0=read_boot(file)/2.0;
  
  // for(size_t iens=0;iens<nens_total;iens++)
  //   cout<<jack_index[0][iens][0]<<endl;

  //cout<<lat_par[0].ml.ave_err()<<endl;

  //alist, zlist
  dbvec_t alist(nbeta),zlist(nbeta);

  //D meson
  dboot_t aMD;
  dboot_t MD_sl_s;

  //output
  dbvec_t output_dM2Pi(ninput_an*nan_syst);
  dbvec_t output_dM2K_QED(ninput_an*nan_syst);
  dbvec_t output_dM2K_QCD_over_minus_two_Deltamud(ninput_an*nan_syst);
  dbvec_t output_dMD_QED(ninput_an*nan_syst);
  dbvec_t output_dMD_QCD(ninput_an*nan_syst);
  dbvec_t output_epsilon(ninput_an*nan_syst);
  dbvec_t output_epsilon_Pi0(ninput_an*nan_syst);
  dbvec_t output_epsilon_K0(ninput_an*nan_syst);

  vector<ave_err_t> v_ave_an_dM2Pi(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QCD_over_minus_two_Deltamud(nan_syst);
  vector<ave_err_t> v_ave_an_dMD_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dMD_QCD(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_Pi0(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_K0(nan_syst);
  
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
	dbvec_t daMD_QED(raw_data.size());
	dbvec_t FVE_daMD(raw_data.size());
	dbvec_t MD(raw_data.size());
	dbvec_t MDs(raw_data.size());
	dbvec_t epsilon_gamma(raw_data.size());
	dbvec_t epsilon_gamma_minusFVE(raw_data.size());
	dbvec_t epsilon_Pi0(raw_data.size());
	dbvec_t epsilon_Pi0_minusFVE(raw_data.size());
	dbvec_t epsilon_K0(raw_data.size());
	dbvec_t epsilon_K0_minusFVE(raw_data.size());
	
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
	    dboot_t aMD_sl_exch=dboot_t(bi,raw_data[iens].D_SL_exch);
	    dboot_t aMD_sl_selftad=dboot_t(bi,raw_data[iens].D_SL_selftad);
	    MD_sl_s=dboot_t(bi,raw_data[iens].D_SL_s);
	    dboot_t MD_sl_p=dboot_t(bi,raw_data[iens].D_SL_p);
	    MD[iens]=aMD/a;
	    
	    daMD_QED[iens]=
	      2*ml[iens]*a*MD_sl_s/Z_QED
	      -(aDeltam_cr_d-aDeltam_cr_u)*MD_sl_p
	      +(sqr(eu)-sqr(ed))*e2*aMD_sl_selftad+(eu-ed)*eu*e2*aMD_sl_exch;
	    FVE_daMD[iens]=FVE_dM(aMD,raw_data[iens].L);
	    
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
	  data_dM2Pi.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						    da2M2Pi[iens]-FVE_da2M2Pi[iens],da2M2Pi[iens]));
	output_dM2Pi[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2Pi(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2Pi,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2Pi_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
        cout<<"                     Epsilon "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						      epsilon_gamma_minusFVE[iens],epsilon_gamma[iens]));
	output_epsilon[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon,lat_par[input_an_id].ml,lat_par[input_an_id].ms,combine("plots/cont_chir_fit_epsilon_gamma_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                      QED K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2K_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						       da2M2K_QED[iens]-FVE_da2M2K[iens],da2M2K_QED[iens]));
	
	output_dM2K_QED[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2K_QED(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QED,lat_par[input_an_id].ml,lat_par[input_an_id].ms,combine("plots/cont_chir_fit_dM2K_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                      QCD K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QCD_over_minus_two_Deltamud;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2K_QCD_over_minus_two_Deltamud.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
									       daM2K_QCD_over_minus_two_Deltamud[iens],daM2K_QCD_over_minus_two_Deltamud[iens]));
	
	output_dM2K_QCD_over_minus_two_Deltamud[input_an_id+an_flag*ninput_an]=cont_chir_fit_dM2K_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QCD_over_minus_two_Deltamud,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2K_QCD_over_minus_two_Deltamud_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         QED D "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dMD_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dMD_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
						      daMD_QED[iens]-FVE_daMD[iens],daMD_QED[iens]));
	
	output_dMD_QED[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_dMD_QED,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dMD_QED_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$(M_{D^+}-M_{D^0})^{QED}[GeV]",1.0,0.0,an_flag,1,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;

	plot_ens_data(combine("plots/MD_an%zu.xmg",input_an_id),ml,MD);
	plot_ens_data(combine("plots/MDs_an%zu.xmg",input_an_id),ml,MDs);
	
	cout<<"                        EpsilonPi0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_Pi0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon_Pi0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							  epsilon_Pi0_minusFVE[iens],epsilon_Pi0[iens]));
	output_epsilon_Pi0[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon_Pi0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_Pi0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_Pi0_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         EpsilonK0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_K0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_epsilon_K0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,raw_data[iens].ibeta,raw_data[iens].L,
							 epsilon_K0_minusFVE[iens],epsilon_K0[iens]));
	output_epsilon_K0[input_an_id+an_flag*ninput_an]=cont_chir_fit_epsilon_K0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_K0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_K0_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag);
	
	cout<<"-----------------------------------------------"<<endl;
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
	dM2K_QCD[input_an_id+an_flag*ninput_an]=dM2K_exp-dboot_t(output_dM2K_QED[input_an_id+an_flag*ninput_an]*1000);
	Deltamud[input_an_id+an_flag*ninput_an]=dboot_t(-dM2K_QCD[input_an_id+an_flag*ninput_an]/output_dM2K_QCD_over_minus_two_Deltamud[input_an_id+an_flag*ninput_an]/2.0);
      }

  /////////////////Deltamud///////////////////
  v_ave_an_Deltamud=ave_analyses(Deltamud);
  cout<<"Deltamud: "<<v_ave_an_Deltamud[7]<<" "<<syst_analysis(v_ave_an_Deltamud)<<endl;

  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t dMD_QCD(raw_data.size());
      
	for(size_t iens=0;iens<raw_data.size();iens++)
	  dMD_QCD[iens]=2*Deltamud[input_an_id+an_flag*ninput_an]*MD_sl_s/1000;

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
      
	output_dMD_QCD[input_an_id+an_flag*ninput_an]=cont_chir_linear_fit(alist,zlist,data_dMD_QCD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dMD_QCD_flag%zu_an%zu.xmg",an_flag,input_an_id),"$$(M_{D^+}-M_{D^0})^{QCD}[GeV]",0.0,1.0,an_flag,0,cov_flag);
      
	cout<<"-----------------------------------------------"<<endl;
      }
  
  ////////////////////output//////////////////
  // v_ave_an_dM2Pi=ave_analyses(output_dM2Pi);
  // cout<<"dM2Pi: "<<v_ave_an_dM2Pi[7]<<" "<<syst_analysis(v_ave_an_dM2Pi)<<endl;
  // v_ave_an_epsilon_Pi0=ave_analyses(output_epsilon_Pi0);
  // cout<<"Epsilon_Pi0: "<<v_ave_an_epsilon_Pi0[7]<<" "<<syst_analysis(v_ave_an_epsilon_Pi0)<<endl;
  // v_ave_an_epsilon=ave_analyses(output_epsilon);
  // cout<<"Epsilon: "<<v_ave_an_epsilon[7]<<" "<<syst_analysis(v_ave_an_epsilon)<<endl;
  v_ave_an_dM2K_QED=ave_analyses(output_dM2K_QED);
  cout<<"dM2K_QED: "<<v_ave_an_dM2K_QED[7]<<" "<<syst_analysis(v_ave_an_dM2K_QED)<<endl;
  // v_ave_an_epsilon_K0=ave_analyses(output_epsilon_K0);
  // cout<<"Epsilon_K0: "<<v_ave_an_epsilon_K0[7]<<" "<<syst_analysis(v_ave_an_epsilon_K0)<<endl;

  
  //////////////////D meson//////////////////
  dbvec_t dMD(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_dMD(nan_syst);

  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      dMD[input_an_id+an_flag*ninput_an]=output_dMD_QCD[input_an_id+an_flag*ninput_an]+output_dMD_QED[input_an_id+an_flag*ninput_an];
  
  v_ave_an_dMD=ave_analyses(dMD);
  cout<<"M_{D^+}-M_{D^0}: "<<v_ave_an_dMD[7]<<endl;
  
  return 0;
}
