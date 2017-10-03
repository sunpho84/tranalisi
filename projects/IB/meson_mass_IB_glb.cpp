#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <meson_mass_IB_fit.hpp>
#include <common.hpp>

vector<ave_err_t> ave_analyses(const dbvec_t &v)
{
  vector<ave_err_t> input_an_ave_err(nan_syst);
  
  for(size_t isyst=0;isyst<nan_syst;isyst++)
    {
      dbvec_t v_an(ninput_an);
      for(size_t inpan=0;inpan<ninput_an;inpan++) v_an[inpan]=v[ind_an({inpan,isyst})];
      input_an_ave_err[isyst]=eq_28_analysis(v_an);
    }
  
  return input_an_ave_err;
}

double syst_analysis(const vector<ave_err_t> &v)
{
  ave_err_t ae;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave();
      ae.ave()+=a;
      ae.err()+=sqr(a);
    }
  ae.ave()/=v.size();
  ae.err()/=v.size();
  ae.err()-=sqr(ae.ave());
  ae.err()=sqrt(fabs(ae.err()));
  
  return ae.err();
}

ave_err_t stat_analysis(const vector<ave_err_t> &v)
{
  ave_err_t ae;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave();
      double e=v[i].err();
      ae.ave()+=a;
      ae.err()+=sqr(e);
    }
  ae.ave()/=v.size();
  ae.err()/=v.size();
  ae.err()=sqrt(fabs(ae.err()));
  
  return ae;
}

void syst_analysis_sep(const vector<ave_err_t> &v)
{
  double db[12];
  
  db[0]=(v[0].ave()-v[1].ave())/2.0;
  db[1]=(v[2].ave()-v[3].ave())/2.0;
  db[2]=(v[4].ave()-v[5].ave())/2.0;
  db[3]=(v[6].ave()-v[7].ave())/2.0;
  db[4]=(v[0].ave()-v[2].ave())/2.0;
  db[5]=(v[1].ave()-v[3].ave())/2.0;
  db[6]=(v[4].ave()-v[6].ave())/2.0;
  db[7]=(v[5].ave()-v[7].ave())/2.0;
  db[8]=(v[0].ave()-v[4].ave())/2.0;
  db[9]=(v[1].ave()-v[5].ave())/2.0;
  db[10]=(v[2].ave()-v[6].ave())/2.0;
  db[11]=(v[3].ave()-v[7].ave())/2.0;
  
  double S2cont,S2chir,S2fse;
  
  S2cont=1.0/24.0*(pow(db[0],2.0)+pow(db[1],2.0)+pow(db[2],2.0)+pow(db[3],2.0)+pow(db[0]+db[1],2.0)+pow(db[0]+db[2],2.0)+pow(db[1]+db[3],2.0)+pow(db[2]+db[3],2.0))+1.0/48.0*(pow(db[0]+db[3],2.0)+pow(db[1]+db[2],2.0));
  
  S2chir=1.0/24.0*(pow(db[8],2.0)+pow(db[9],2.0)+pow(db[10],2.0)+pow(db[11],2.0)+pow(db[8]+db[9],2.0)+pow(db[8]+db[10],2.0)+pow(db[9]+db[11],2.0)+pow(db[10]+db[11],2.0))+1.0/48.0*(pow(db[8]+db[11],2.0)+pow(db[9]+db[10],2.0));
  
  S2fse=1.0/24.0*(pow(db[4],2.0)+pow(db[5],2.0)+pow(db[6],2.0)+pow(db[7],2.0)+pow(db[4]+db[5],2.0)+pow(db[4]+db[6],2.0)+pow(db[5]+db[7],2.0)+pow(db[6]+db[7],2.0))+1.0/48.0*(pow(db[4]+db[7],2.0)+pow(db[5]+db[6],2.0));

  cout<<"cont: "<<sqrt(S2cont)<<endl;
  cout<<"chir: "<<sqrt(S2chir)<<endl;
  cout<<"fse: "<<sqrt(S2fse)<<endl;
}

class ens_data_t
{
public:
  size_t iult,ibeta,L,useforL,useforL1,useforL2;
  double aml,ams,amc;
  string path;
  
  djack_t pi_mass,pi_SL_exch,pi_SL_selftad,pi_SL_s,pi_SL_p,k_mass,k_SL_exch,k_SL_selftad,k_SL_selftad_revins,k_SL_s,k_A_s,k_SL_s_revins,k_SL_p,k_SL_p_revins,D_mass,D_SL_exch,D_SL_selftad,D_SL_selftad_revins,D_SL_s,D_SL_s_revins,D_SL_p,D_SL_p_revins,Ds_mass,Ds_SL_exch,Ds_SL_selftad,Ds_SL_selftad_revins,Ds_SL_s,Ds_SL_s_revins,Ds_SL_p,Ds_SL_p_revins;
  djack_t deltam_cr;
};

vector<ens_data_t> raw_data;

//! plot the three ensemble separately
template <class Tx,class Ty> void plot_ens_data(grace_file_t &file,const Tx &x,const Ty &y)
{
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    {
      file.new_data_set();
      for(size_t iens=0;iens<raw_data.size();iens++)
	if(raw_data[iens].ibeta==ibeta)
	  file.write_ave_err(x[iens].ave(),y[iens].ave_err());
    }
}

template <class Tx,class Ty> void plot_ens_data(string path,const Tx &x,const Ty &y)
{
  grace_file_t file(path);
  plot_ens_data(file,x,y);
}

int main(int narg,char **arg)
{
  ind_an.set_ranges({{"Input",ninput_an},{"Syst",nan_syst}});
  
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
      temp.useforL1=input.read<size_t>("useforL1");
      temp.useforL2=input.read<size_t>("useforL2");
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
      obs_k_file.bin_read(temp.k_SL_s_revins);
      obs_k_file.bin_read(temp.k_A_s);

      //read the observable of the D meson
      raw_file_t obs_D_file(combine("%s/D_obs",temp.path.c_str()),"r");
      obs_D_file.bin_read(temp.D_mass);
      obs_D_file.bin_read(temp.D_SL_exch);
      obs_D_file.bin_read(temp.D_SL_selftad);
      obs_D_file.bin_read(temp.D_SL_s);
      obs_D_file.bin_read(temp.D_SL_p);
      obs_D_file.bin_read(temp.D_SL_selftad_revins);
      obs_D_file.bin_read(temp.D_SL_p_revins);
      obs_D_file.bin_read(temp.D_SL_s_revins);

      //read the observable of the Ds meson
      raw_file_t obs_Ds_file(combine("%s/Ds_obs",temp.path.c_str()),"r");
      obs_Ds_file.bin_read(temp.Ds_mass);
      obs_Ds_file.bin_read(temp.Ds_SL_exch);
      obs_Ds_file.bin_read(temp.Ds_SL_selftad);
      obs_Ds_file.bin_read(temp.Ds_SL_s);
      obs_Ds_file.bin_read(temp.Ds_SL_p);
      obs_Ds_file.bin_read(temp.Ds_SL_selftad_revins);
      obs_Ds_file.bin_read(temp.Ds_SL_p_revins);
      obs_Ds_file.bin_read(temp.Ds_SL_s_revins);

      //read deltam_cr (ud)
      raw_file_t deltam_cr_file(combine("%s/ud_fit_deltam_cr",temp.path.c_str()),"r");
      deltam_cr_file.bin_read(temp.deltam_cr);

      //store in the raw_data vector
      raw_data.push_back(temp);
    }

  //alist, zlist
  dbvec_t alist(nbeta),zlist(nbeta);
  const vector<string> beta_list={"1.90","1.95","2.10"};
  
  //D meson + scalar densities
  dbvec_t MPi_sl_s(ninput_an*raw_data.size());
  dbvec_t MK_sl_s(ninput_an*raw_data.size());
  dbvec_t DG2K_fr_G2K(ninput_an*raw_data.size());
  dbvec_t MK_sl_s_revins(ninput_an*raw_data.size());
  dbvec_t aDeltam_cr_u(ninput_an*raw_data.size());
  dbvec_t aDeltam_cr_d(ninput_an*raw_data.size());
  dbvec_t MPi(ninput_an*raw_data.size());
  dbvec_t MK(ninput_an*raw_data.size());
  dbvec_t MD(ninput_an*raw_data.size());
  dbvec_t MDs(ninput_an*raw_data.size());
  dbvec_t aMD(ninput_an*raw_data.size());
  dbvec_t aMDs(ninput_an*raw_data.size());
  dbvec_t aMD_sl_exch(ninput_an*raw_data.size());
  dbvec_t aMD_sl_selftad(ninput_an*raw_data.size());
  dbvec_t MD_sl_s(ninput_an*raw_data.size());
  dbvec_t MD_sl_s_revins(ninput_an*raw_data.size());
  dbvec_t MD_sl_p(ninput_an*raw_data.size());
  dbvec_t aMD_sl_selftad_revins(ninput_an*raw_data.size());
  dbvec_t MD_sl_p_revins(ninput_an*raw_data.size());
  dbvec_t MDs_sl_s_revins(ninput_an*raw_data.size());
  dbvec_t Z_QED(ninput_an*raw_data.size());
  
  //output
  dbvec_t output_dM2Pi(ninput_an*nan_syst);
  dbvec_t output_dM2K_QED(ninput_an*nan_syst);
  dbvec_t output_dM2K_QCD_over_minus_two_Deltamud(ninput_an*nan_syst);
  dbvec_t output_dM2D_QED(ninput_an*nan_syst);
  dbvec_t output_dM2D_QCD(ninput_an*nan_syst);
  dbvec_t output_dM2D_QCD_ind(ninput_an*nan_syst);
  dbvec_t output_sMD(ninput_an*nan_syst);
  dbvec_t output_MDs(ninput_an*nan_syst);
  dbvec_t output_epsilon(ninput_an*nan_syst);
  dbvec_t output_epsilon_Pi0(ninput_an*nan_syst);
  dbvec_t output_epsilon_K0(ninput_an*nan_syst);
  dbvec_t output_M2Pi0g(ninput_an*nan_syst);
  dbvec_t output_M2K0g(ninput_an*nan_syst);
  dbvec_t output_M2uug(ninput_an*nan_syst);
  dbvec_t output_M2ddg(ninput_an*nan_syst);
  dbvec_t output_DeltaFK_over_Deltamud(ninput_an*nan_syst);

  vector<ave_err_t> v_ave_an_dM2Pi(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QCD_over_minus_two_Deltamud(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_QED(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_QCD(nan_syst);
  vector<ave_err_t> v_ave_an_dM2D_QCD_ind(nan_syst);
  vector<ave_err_t> v_ave_an_sMD(nan_syst);
  vector<ave_err_t> v_ave_an_MDs(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_Pi0(nan_syst);
  vector<ave_err_t> v_ave_an_epsilon_K0(nan_syst);
  vector<ave_err_t> v_ave_an_M2Pi0g(nan_syst);
  vector<ave_err_t> v_ave_an_M2K0g(nan_syst);
  vector<ave_err_t> v_ave_an_M2uug(nan_syst);
  vector<ave_err_t> v_ave_an_M2ddg(nan_syst);
  vector<ave_err_t> v_ave_an_DeltaFK_over_Deltamud(nan_syst);
 
  bool cov_flag=false;

  //loop over analysis flags and input scale determination
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t ml(raw_data.size());
	dbvec_t ms(raw_data.size());
	dbvec_t a2M2Pi(raw_data.size());
	dbvec_t da2M2Pi(raw_data.size());
	dbvec_t FVE_da2M2Pi(raw_data.size());
	dbvec_t da2M2K_QED(raw_data.size());
	dbvec_t FVE_da2M2K(raw_data.size());
	dbvec_t daM2K_QCD_over_minus_two_Deltamud(raw_data.size());
	dbvec_t da2M2D_QED(raw_data.size());
	dbvec_t FVE_da2M2D(raw_data.size());
	dbvec_t epsilon_gamma(raw_data.size());
	dbvec_t epsilon_gamma_minusFVE(raw_data.size());
	dbvec_t epsilon_Pi0(raw_data.size());
	dbvec_t epsilon_Pi0_minusFVE(raw_data.size());
	dbvec_t epsilon_K0(raw_data.size());
	dbvec_t epsilon_K0_minusFVE(raw_data.size());
	dbvec_t a2M2Pi0g(raw_data.size());
	dbvec_t a2M2K0g(raw_data.size());
	dbvec_t a2M2uug(raw_data.size());
	dbvec_t a2M2ddg(raw_data.size());
	dbvec_t DeltaFK_over_Deltamud(raw_data.size());
	
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    size_t ens_id=raw_data[iens].iult;
	    size_t ibeta=raw_data[iens].ibeta;
	    
	    dboot_t a=1.0/lat_par[input_an_id].ainv[ibeta];
	    dboot_t Lphys=raw_data[iens].L*a;
	    ml[iens]=raw_data[iens].aml/lat_par[input_an_id].Z[ibeta]/a;
	    ms[iens]=raw_data[iens].ams/lat_par[input_an_id].Z[ibeta]/a;
	    
	    boot_init_t &bi=jack_index[input_an_id][ens_id];
	    
	    Z_QED[iens+input_an_id*raw_data.size()]=1.0/(e2*lat_par[input_an_id].Z[ibeta]*(6.0*log(mu_MS*a)-22.596)/(16.0*sqr(M_PI)));
	    aDeltam_cr_u[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(eu);
	    aDeltam_cr_d[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(ed);
	    
	    dboot_t aMPi=dboot_t(bi,raw_data[iens].pi_mass);
	    dboot_t aMPi_sl_exch=dboot_t(bi,raw_data[iens].pi_SL_exch);
	    dboot_t aMPi_sl_selftad=dboot_t(bi,raw_data[iens].pi_SL_selftad);
	    MPi_sl_s[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].pi_SL_s);
	    dboot_t MPi_sl_p=dboot_t(bi,raw_data[iens].pi_SL_p);
	    a2M2Pi[iens]=aMPi*aMPi;
	    da2M2Pi[iens]=aMPi*sqr(eu-ed)*e2*aMPi_sl_exch;
	    FVE_da2M2Pi[iens]=FVE_M2(aMPi,raw_data[iens].L);
	    MPi[iens+input_an_id*raw_data.size()]=aMPi/a;
	    
	    dboot_t aMK=dboot_t(bi,raw_data[iens].k_mass);
	    dboot_t aMK_sl_exch=dboot_t(bi,raw_data[iens].k_SL_exch);
	    dboot_t aMK_sl_selftad=dboot_t(bi,raw_data[iens].k_SL_selftad);
	    MK_sl_s[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].k_SL_s);
	    MK_sl_s_revins[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].k_SL_s_revins);
	    dboot_t MK_sl_p=dboot_t(bi,raw_data[iens].k_SL_p);
	    dboot_t aMK_sl_selftad_revins=dboot_t(bi,raw_data[iens].k_SL_selftad_revins);
	    dboot_t MK_sl_p_revins=dboot_t(bi,raw_data[iens].k_SL_p_revins);
	    DG2K_fr_G2K[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].k_A_s);
	    
	    dboot_t daMK_QED=
	      -(aDeltam_cr_u[iens+input_an_id*raw_data.size()]-aDeltam_cr_d[iens+input_an_id*raw_data.size()])*MK_sl_p-2.0*ml[iens]*a*MK_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]*2.0/(sqr(ed)-sqr(eu)))
	      +(sqr(eu)-sqr(ed))*e2*(aMK_sl_exch-aMK_sl_selftad);
	    da2M2K_QED[iens]=daMK_QED*2.0*aMK;
	    FVE_da2M2K[iens]=FVE_M2(aMK,raw_data[iens].L);
	    daM2K_QCD_over_minus_two_Deltamud[iens]=2.0*aMK*MK_sl_s[iens+input_an_id*raw_data.size()];
	    MK[iens+input_an_id*raw_data.size()]=aMK/a;
	    
	    aMD[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_mass);
	    aMD_sl_exch[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_exch);
	    aMD_sl_selftad[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_selftad);
	    MD_sl_s[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_s);
	    MD_sl_s_revins[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_s_revins);
	    MD_sl_p[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_p);
	    aMD_sl_selftad_revins[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_selftad_revins);
	    MD_sl_p_revins[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].D_SL_p_revins);
	    MD[iens+input_an_id*raw_data.size()]=aMD[iens+input_an_id*raw_data.size()]/a;
	    
	    dboot_t daMD_QED=
	      2.0*ml[iens]*a*MD_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]*2.0/(sqr(ed)-sqr(eu)))
	      -(aDeltam_cr_d[iens+input_an_id*raw_data.size()]-aDeltam_cr_u[iens+input_an_id*raw_data.size()])*MD_sl_p[iens+input_an_id*raw_data.size()]
	      +(sqr(eu)-sqr(ed))*e2*aMD_sl_selftad[iens+input_an_id*raw_data.size()]+(eu-ed)*eu*e2*aMD_sl_exch[iens+input_an_id*raw_data.size()];
	    da2M2D_QED[iens]=daMD_QED*2.0*aMD[iens+input_an_id*raw_data.size()];
	    FVE_da2M2D[iens]=FVE_M2(aMD[iens+input_an_id*raw_data.size()],raw_data[iens].L);
	    
	    aMDs[iens+input_an_id*raw_data.size()]=dboot_t(bi,raw_data[iens].Ds_mass);
	    MDs[iens+input_an_id*raw_data.size()]=aMDs[iens+input_an_id*raw_data.size()]/a;
	    
	    epsilon_gamma[iens]=(da2M2K_QED[iens]/da2M2Pi[iens])-1.0;
	    epsilon_gamma_minusFVE[iens]=((da2M2K_QED[iens]-FVE_da2M2K[iens])/(da2M2Pi[iens]-FVE_da2M2Pi[iens]))-1.0;
	    
	    dboot_t num_epsilon_Pi0=2.0*aMPi*(-(sqr(eu)+sqr(ed))*e2*(aMPi_sl_exch/2.0+aMPi_sl_selftad)-(aDeltam_cr_u[iens+input_an_id*raw_data.size()]+aDeltam_cr_d[iens+input_an_id*raw_data.size()])*MPi_sl_p+2.0*ml[iens]*a*MPi_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]*2.0/(sqr(ed)+sqr(eu))));
	    epsilon_Pi0[iens]=num_epsilon_Pi0/da2M2Pi[iens];
	    epsilon_Pi0_minusFVE[iens]=num_epsilon_Pi0/(da2M2Pi[iens]-FVE_da2M2Pi[iens]);
	    
	    dboot_t num_epsilon_K0=2.0*aMK*(-sqr(ed)*e2*(aMK_sl_exch+aMK_sl_selftad_revins+aMK_sl_selftad)-aDeltam_cr_d[iens+input_an_id*raw_data.size()]*(MK_sl_p+MK_sl_p_revins)+ml[iens]*a*MK_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(ed))+ms[iens]*a*MK_sl_s_revins[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(ed)));
	    epsilon_K0[iens]=num_epsilon_K0/da2M2Pi[iens];
	    epsilon_K0_minusFVE[iens]=num_epsilon_K0/(da2M2Pi[iens]-FVE_da2M2Pi[iens]);
	    a2M2Pi0g[iens]=num_epsilon_Pi0;
	    a2M2K0g[iens]=num_epsilon_K0;

	    //////////uu; dd mesons/////////
	    dboot_t uu_meson=2.0*aMPi*(-sqr(eu)*e2*(aMPi_sl_exch+2.0*aMPi_sl_selftad)-2.0*aDeltam_cr_u[iens+input_an_id*raw_data.size()]*MPi_sl_p+2.0*ml[iens]*a*MPi_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(eu)));
	    dboot_t dd_meson=2.0*aMPi*(-sqr(ed)*e2*(aMPi_sl_exch+2.0*aMPi_sl_selftad)-2.0*aDeltam_cr_d[iens+input_an_id*raw_data.size()]*MPi_sl_p+2.0*ml[iens]*a*MPi_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(ed)));	    
	    a2M2uug[iens]=uu_meson;
	    a2M2ddg[iens]=dd_meson;

	    //////////DeltaFK////////////
	    DeltaFK_over_Deltamud[iens]=-1.0/(raw_data[iens].ams+raw_data[iens].aml)+DG2K_fr_G2K[iens+input_an_id*raw_data.size()]/2.0+2.0*MK_sl_s[iens+input_an_id*raw_data.size()]/aMK[iens+input_an_id*raw_data.size()];
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
	/*
	cout<<"                         Pi "<<endl;
	cout<<endl;
	
	//data to fit
	vector<cont_chir_fit_data_t> data_dM2Pi;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_dM2Pi.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						      da2M2Pi[iens]-FVE_da2M2Pi[iens],da2M2Pi[iens]));
	
	output_dM2Pi[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2Pi(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2Pi,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2Pi_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
		
	cout<<"                      QED K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(FSE_an(an_flag)==1 or raw_data[iens].useforL)
	    data_dM2K_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							 da2M2K_QED[iens]-FVE_da2M2K[iens],da2M2K_QED[iens]));
	
	output_dM2K_QED[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2K_QED(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QED,lat_par[input_an_id].ml,lat_par[input_an_id].ms,combine("plots/cont_chir_fit_dM2K_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                      QCD K "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2K_QCD_over_minus_two_Deltamud;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2K_QCD_over_minus_two_Deltamud.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
									       daM2K_QCD_over_minus_two_Deltamud[iens],daM2K_QCD_over_minus_two_Deltamud[iens]));
	
	output_dM2K_QCD_over_minus_two_Deltamud[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2K_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2K_QCD_over_minus_two_Deltamud,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2K_QCD_over_minus_two_Deltamud_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         M2Pi0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2Pi0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_M2Pi0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						       a2M2Pi0g[iens],a2M2Pi0g[iens]));
	
	output_M2Pi0g[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2Pi0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2Pi0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2Pi0g_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         M2uug "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2uug;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_M2uug.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						       a2M2uug[iens],a2M2uug[iens]));
	
	output_M2uug[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2Pi0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2uug,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2uug_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         M2ddg "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2ddg;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_M2ddg.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						       a2M2ddg[iens],a2M2ddg[iens]));
	
	output_M2ddg[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2Pi0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2ddg,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2ddg_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;

	cout<<"                         M2K0g "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_M2K0g;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_M2K0g.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							 a2M2K0g[iens],a2M2K0g[iens]));
	
	output_M2K0g[ind_an({input_an_id,an_flag})]=cont_chir_fit_M2K0g(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_M2K0g,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_M2K0g_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                     Epsilon "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_epsilon.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						      epsilon_gamma_minusFVE[iens],epsilon_gamma[iens]));
	
	output_epsilon[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_gamma_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                        EpsilonPi0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_Pi0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_epsilon_Pi0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							  epsilon_Pi0_minusFVE[iens],epsilon_Pi0[iens]));
	
	output_epsilon_Pi0[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon_Pi0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_Pi0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_Pi0_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	cout<<"                         EpsilonK0 "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_epsilon_K0;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_epsilon_K0.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							 epsilon_K0_minusFVE[iens],epsilon_K0[iens]));
	
	output_epsilon_K0[ind_an({input_an_id,an_flag})]=cont_chir_fit_epsilon_K0(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_epsilon_K0,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_epsilon_K0_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         QED D "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_dM2D_QED;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_dM2D_QED.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						      da2M2D_QED[iens]-FVE_da2M2D[iens],da2M2D_QED[iens]));
	
	output_dM2D_QED[ind_an({input_an_id,an_flag})]=cont_chir_linear_fit(alist,zlist,data_dM2D_QED,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QED_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),"$$[M^2_{D^+}-M^2_{D^0})_{QED} [GeV^2]",2.0,0.0,an_flag,1,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	*/
	cout<<"                      DeltaFK "<<endl;
	cout<<endl;
	
	vector<cont_chir_fit_data_t> data_DeltaFK_over_Deltamud;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  data_DeltaFK_over_Deltamud.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
									       DeltaFK_over_Deltamud[iens],DeltaFK_over_Deltamud[iens]));
	if(input_an_id==7) fit_debug=1;	
	output_DeltaFK_over_Deltamud[ind_an({input_an_id,an_flag})]=cont_chir_fit_DeltaFK(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_DeltaFK_over_Deltamud,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_DeltaFK_over_Deltamud_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	cout<<endl;
	
	dbvec_t MPi_o(raw_data.size());
	dbvec_t MK_o(raw_data.size());
	dbvec_t MD_o(raw_data.size());
	dbvec_t MDs_o(raw_data.size());
	dbvec_t MPi_sl_s_o(raw_data.size());
	dbvec_t MK_sl_s_o(raw_data.size());
	dbvec_t MK_sl_s_revins_o(raw_data.size());
	dbvec_t MD_sl_s_o(raw_data.size());
	dbvec_t MD_sl_s_revins_o(raw_data.size());
	    
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    MPi_o[iens]=MPi[iens+input_an_id*raw_data.size()];
	    MK_o[iens]=MK[iens+input_an_id*raw_data.size()];
	    MD_o[iens]=MD[iens+input_an_id*raw_data.size()];
	    MDs_o[iens]=MDs[iens+input_an_id*raw_data.size()];
	    MPi_sl_s_o[iens]=MPi_sl_s[iens+input_an_id*raw_data.size()];
	    MK_sl_s_o[iens]=MK_sl_s[iens+input_an_id*raw_data.size()];
	    MK_sl_s_revins_o[iens]=MK_sl_s_revins[iens+input_an_id*raw_data.size()];
	    MD_sl_s_o[iens]=MD_sl_s[iens+input_an_id*raw_data.size()];
	    MD_sl_s_revins_o[iens]=MD_sl_s_revins[iens+input_an_id*raw_data.size()];
	  }

	plot_ens_data(combine("plots/MPi_an%zu.xmg",input_an_id),ml,MPi_o);
	plot_ens_data(combine("plots/MK_an%zu.xmg",input_an_id),ml,MK_o);
	plot_ens_data(combine("plots/MD_an%zu.xmg",input_an_id),ml,MD_o);
	plot_ens_data(combine("plots/MDs_an%zu.xmg",input_an_id),ml,MDs_o);
	plot_ens_data(combine("plots/MPi_sl_s_an%zu.xmg",input_an_id),ml,MPi_sl_s_o);
	plot_ens_data(combine("plots/MK_sl_s_an%zu.xmg",input_an_id),ml,MK_sl_s_o);
	plot_ens_data(combine("plots/MK_sl_s_revins_an%zu.xmg",input_an_id),ml,MK_sl_s_revins_o);
	plot_ens_data(combine("plots/MD_sl_s_an%zu.xmg",input_an_id),ml,MD_sl_s_o);
	plot_ens_data(combine("plots/MD_sl_s_revins_an%zu.xmg",input_an_id),ml,MD_sl_s_revins_o);
		
      }
  /*
  ///////////////////masses//////////////////
  for(size_t iens=0;iens<raw_data.size();iens++)
    {
      dbvec_t v_ens_MPi(ninput_an);
      dbvec_t v_ens_MK(ninput_an);
      dbvec_t v_ens_MD(ninput_an);
      dbvec_t v_ens_MDs(ninput_an);

      for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
	{
	  v_ens_MPi[input_an_id]=MPi[iens+input_an_id*raw_data.size()];
	  v_ens_MK[input_an_id]=MK[iens+input_an_id*raw_data.size()];
	  v_ens_MD[input_an_id]=MD[iens+input_an_id*raw_data.size()];
	  v_ens_MDs[input_an_id]=MDs[iens+input_an_id*raw_data.size()];	  
	}

      cout<<"iens: "<<iens<<endl;
      cout<<"MPi: "<<eq_28_analysis(v_ens_MPi)<<endl;
      cout<<"MK: "<<eq_28_analysis(v_ens_MK)<<endl;
      cout<<"MD: "<<eq_28_analysis(v_ens_MD)<<endl;
      cout<<"MDs: "<<eq_28_analysis(v_ens_MDs)<<endl;
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

  v_ave_an_M2uug=ave_analyses(output_M2uug);
  cout<<"M2uug: "<<stat_analysis(v_ave_an_M2uug)<<" "<<syst_analysis(v_ave_an_M2uug)<<endl;
  syst_analysis_sep(v_ave_an_M2uug);
  for(size_t i=0;i<8;i++)
    cout<<"an_M2uug: "<<v_ave_an_M2uug[i]<<endl;

  v_ave_an_M2ddg=ave_analyses(output_M2ddg);
  cout<<"M2ddg: "<<stat_analysis(v_ave_an_M2ddg)<<" "<<syst_analysis(v_ave_an_M2ddg)<<endl;
  syst_analysis_sep(v_ave_an_M2ddg);
  for(size_t i=0;i<8;i++)
    cout<<"an_M2ddg: "<<v_ave_an_M2ddg[i]<<endl;
  
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
  */
  v_ave_an_DeltaFK_over_Deltamud=ave_analyses(output_DeltaFK_over_Deltamud);
  cout<<"DeltaFK_over_Deltamud: "<<stat_analysis(v_ave_an_DeltaFK_over_Deltamud)<<" "<<syst_analysis(v_ave_an_DeltaFK_over_Deltamud)<<endl;
  syst_analysis_sep(v_ave_an_DeltaFK_over_Deltamud);
  for(size_t i=0;i<8;i++)
    cout<<"an_DeltaFK_over_Deltamud: "<<v_ave_an_DeltaFK_over_Deltamud[i]<<endl;
  /*
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
	epsilon_Pi0_ind[ind_an({input_an_id,an_flag})]=dboot_t(output_M2Pi0g[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]);
	epsilon_K0_ind[ind_an({input_an_id,an_flag})]=dboot_t(output_M2K0g[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]);
	epsilon_gamma_ind[ind_an({input_an_id,an_flag})]=dboot_t(output_dM2K_QED[ind_an({input_an_id,an_flag})]/output_dM2Pi[ind_an({input_an_id,an_flag})]-1.0);
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

  double M_pion_QCD=mpi0-v_ave_an_epsilon_Pi0[7].ave()*(mpip-mpi0);
  double M_pion_QCD_ind=mpi0-v_ave_an_epsilon_Pi0_ind[7].ave()*(mpip-mpi0);

  cout<<"pion QCD: "<<M_pion_QCD<<endl;
  cout<<"pion QCD_ind: "<<M_pion_QCD_ind<<endl;

  ///////////////////Deltamud + Q2 Colangelo/////////////////
  dbvec_t dM2K_QED_ind(ninput_an*nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QED_ind[ind_an({input_an_id,an_flag})]=dboot_t((output_epsilon[ind_an({input_an_id,an_flag})]+1.0)*output_dM2Pi[ind_an({input_an_id,an_flag})]);
      }
  
  dbvec_t dM2K_QCD(ninput_an*nan_syst);
  dbvec_t dM2K_QCD_ind(ninput_an*nan_syst);
  dbvec_t Deltamud(ninput_an*nan_syst);
  dbvec_t Deltamud_ind(ninput_an*nan_syst);
  dbvec_t Q2_col(ninput_an*nan_syst);
  dbvec_t Q2_col_ind(ninput_an*nan_syst);
  dbvec_t DeltaFK(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QCD(nan_syst);
  vector<ave_err_t> v_ave_an_dM2K_QCD_ind(nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud(nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud_ind(nan_syst);
  vector<ave_err_t> v_ave_an_Q2_col(nan_syst);
  vector<ave_err_t> v_ave_an_Q2_col_ind(nan_syst);
  vector<ave_err_t> v_ave_an_DeltaFK(nan_syst);
  dboot_t dM2K_exp,M2Pi_exp,M2K_exp;
  for(size_t iboot=0;iboot<nboots;iboot++)
    {
      dM2K_exp[iboot]=-3.903;
      M2Pi_exp[iboot]=pow(134.98,2.0);
      M2K_exp[iboot]=pow(494.2,2.0);
    }
  
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dM2K_QCD[ind_an({input_an_id,an_flag})]=dM2K_exp-dboot_t(output_dM2K_QED[ind_an({input_an_id,an_flag})]*1000.0);
	dM2K_QCD_ind[ind_an({input_an_id,an_flag})]=dM2K_exp-dboot_t(dM2K_QED_ind[ind_an({input_an_id,an_flag})]*1000.0);
	Deltamud[ind_an({input_an_id,an_flag})]=dboot_t(-dM2K_QCD[ind_an({input_an_id,an_flag})]/output_dM2K_QCD_over_minus_two_Deltamud[ind_an({input_an_id,an_flag})]/2.0);
	Deltamud_ind[ind_an({input_an_id,an_flag})]=dboot_t(-dM2K_QCD_ind[ind_an({input_an_id,an_flag})]/output_dM2K_QCD_over_minus_two_Deltamud[ind_an({input_an_id,an_flag})]/2.0);
	Q2_col[ind_an({input_an_id,an_flag})]=dboot_t(-M2K_exp*(M2K_exp-M2Pi_exp)/M2Pi_exp/(dM2K_QCD[ind_an({input_an_id,an_flag})]*1000.0));
	Q2_col_ind[ind_an({input_an_id,an_flag})]=dboot_t(-M2K_exp*(M2K_exp-M2Pi_exp)/M2Pi_exp/(dM2K_QCD_ind[ind_an({input_an_id,an_flag})]*1000.0));
	DeltaFK[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud[ind_an({input_an_id,an_flag})]*output_DeltaFK_over_Deltamud[ind_an({input_an_id,an_flag})]);
      }

  v_ave_an_dM2K_QCD=ave_analyses(dM2K_QCD);
  cout<<"dM2K_QCD: "<<stat_analysis(v_ave_an_dM2K_QCD)<<" "<<syst_analysis(v_ave_an_dM2K_QCD)<<endl;
  syst_analysis_sep(v_ave_an_dM2K_QCD);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QCD: "<<v_ave_an_dM2K_QCD[i]<<endl;

  v_ave_an_dM2K_QCD_ind=ave_analyses(dM2K_QCD_ind);
  cout<<"dM2K_QCD_ind: "<<stat_analysis(v_ave_an_dM2K_QCD_ind)<<" "<<syst_analysis(v_ave_an_dM2K_QCD_ind)<<endl;
  syst_analysis_sep(v_ave_an_dM2K_QCD_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_dM2K_QCD_ind: "<<v_ave_an_dM2K_QCD_ind[i]<<endl;
    
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
  
  v_ave_an_Q2_col=ave_analyses(Q2_col);
  cout<<"Q2_col: "<<stat_analysis(v_ave_an_Q2_col)<<" "<<syst_analysis(v_ave_an_Q2_col)<<endl;
  syst_analysis_sep(v_ave_an_Q2_col);
  for(size_t i=0;i<8;i++)
    cout<<"an_Q2_col: "<<v_ave_an_Q2_col[i]<<endl;

  v_ave_an_Q2_col_ind=ave_analyses(Q2_col_ind);
  cout<<"Q2_col_ind: "<<stat_analysis(v_ave_an_Q2_col_ind)<<" "<<syst_analysis(v_ave_an_Q2_col_ind)<<endl;
  syst_analysis_sep(v_ave_an_Q2_col_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Q2_col_ind: "<<v_ave_an_Q2_col_ind[i]<<endl;

  v_ave_an_DeltaFK=ave_analyses(DeltaFK);
  cout<<"DeltaFK: "<<stat_analysis(v_ave_an_DeltaFK)<<" "<<syst_analysis(v_ave_an_DeltaFK)<<endl;
  syst_analysis_sep(v_ave_an_DeltaFK);
  for(size_t i=0;i<8;i++)
    cout<<"an_DeltaFK: "<<v_ave_an_DeltaFK[i]<<endl;
  
  /////////////////R,Q2,mu/md/////////////////////
  dbvec_t R(ninput_an*nan_syst);
  dbvec_t R_ind(ninput_an*nan_syst);
  dbvec_t Q2(ninput_an*nan_syst);
  dbvec_t Q2_ind(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud_ind(ninput_an*nan_syst);
  dbvec_t ratio_mu_md(ninput_an*nan_syst);
  dbvec_t ratio_mu_md_ind(ninput_an*nan_syst);
  dbvec_t Deltamud_col(ninput_an*nan_syst);
  dbvec_t Deltamud_col_ind(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud_col(ninput_an*nan_syst);
  dbvec_t Deltamud_over_mud_col_ind(ninput_an*nan_syst);
  dbvec_t ratio_mu_md_col(ninput_an*nan_syst);
  dbvec_t ratio_mu_md_col_ind(ninput_an*nan_syst);
  dbvec_t mu_q(ninput_an*nan_syst);
  dbvec_t md_q(ninput_an*nan_syst);
  vector<ave_err_t> v_ave_an_R(nan_syst);
  vector<ave_err_t> v_ave_an_R_ind(nan_syst);
  vector<ave_err_t> v_ave_an_Q2(nan_syst);
  vector<ave_err_t> v_ave_an_Q2_ind(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md_ind(nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud_col(nan_syst);
  vector<ave_err_t> v_ave_an_Deltamud_col_ind(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md_col(nan_syst);
  vector<ave_err_t> v_ave_an_ratio_mu_md_col_ind(nan_syst);
  vector<ave_err_t> v_ave_an_mu(nan_syst);
  vector<ave_err_t> v_ave_an_md(nan_syst);
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	R[ind_an({input_an_id,an_flag})]=dboot_t((lat_par[input_an_id].ms-lat_par[input_an_id].ml)/(2.0*Deltamud[ind_an({input_an_id,an_flag})]*1.0e-3));
	R_ind[ind_an({input_an_id,an_flag})]=dboot_t((lat_par[input_an_id].ms-lat_par[input_an_id].ml)/(2.0*Deltamud_ind[ind_an({input_an_id,an_flag})]*1.0e-3));
	Q2[ind_an({input_an_id,an_flag})]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Deltamud[ind_an({input_an_id,an_flag})]*lat_par[input_an_id].ml*1.0e-3));
	Q2_ind[ind_an({input_an_id,an_flag})]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Deltamud_ind[ind_an({input_an_id,an_flag})]*lat_par[input_an_id].ml*1.0e-3));
	Deltamud_over_mud[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud[ind_an({input_an_id,an_flag})]*1.0e-3/lat_par[input_an_id].ml);
	Deltamud_over_mud_ind[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud_ind[ind_an({input_an_id,an_flag})]*1.0e-3/lat_par[input_an_id].ml);
	ratio_mu_md[ind_an({input_an_id,an_flag})]=dboot_t((1.0-Deltamud_over_mud[ind_an({input_an_id,an_flag})])/(1.0+Deltamud_over_mud[ind_an({input_an_id,an_flag})]));
	ratio_mu_md_ind[ind_an({input_an_id,an_flag})]=dboot_t((1.0-Deltamud_over_mud_ind[ind_an({input_an_id,an_flag})])/(1.0+Deltamud_over_mud_ind[ind_an({input_an_id,an_flag})]));
	Deltamud_col[ind_an({input_an_id,an_flag})]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Q2_col[ind_an({input_an_id,an_flag})]*lat_par[input_an_id].ml*1.0e-3));
	Deltamud_col_ind[ind_an({input_an_id,an_flag})]=dboot_t((pow(lat_par[input_an_id].ms,2.0)-pow(lat_par[input_an_id].ml,2.0))/(4.0*Q2_col_ind[ind_an({input_an_id,an_flag})]*lat_par[input_an_id].ml*1.0e-3));
	Deltamud_over_mud_col[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud_col[ind_an({input_an_id,an_flag})]*1.0e-3/lat_par[input_an_id].ml);
	Deltamud_over_mud_col_ind[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud_col_ind[ind_an({input_an_id,an_flag})]*1.0e-3/lat_par[input_an_id].ml);
	ratio_mu_md_col[ind_an({input_an_id,an_flag})]=dboot_t((1.0-Deltamud_over_mud_col[ind_an({input_an_id,an_flag})])/(1.0+Deltamud_over_mud_col[ind_an({input_an_id,an_flag})]));
	ratio_mu_md_col_ind[ind_an({input_an_id,an_flag})]=dboot_t((1.0-Deltamud_over_mud_col_ind[ind_an({input_an_id,an_flag})])/(1.0+Deltamud_over_mud_col_ind[ind_an({input_an_id,an_flag})]));
	mu_q[ind_an({input_an_id,an_flag})]=dboot_t(ratio_mu_md[ind_an({input_an_id,an_flag})]*Deltamud[ind_an({input_an_id,an_flag})]*2.0/(1.0-ratio_mu_md[ind_an({input_an_id,an_flag})]));
	md_q[ind_an({input_an_id,an_flag})]=dboot_t(Deltamud[ind_an({input_an_id,an_flag})]*2.0/(1.0-ratio_mu_md[ind_an({input_an_id,an_flag})]));
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

  v_ave_an_Deltamud_col=ave_analyses(Deltamud_col);
  cout<<"Deltamud_col: "<<stat_analysis(v_ave_an_Deltamud_col)<<" "<<syst_analysis(v_ave_an_Deltamud_col)<<endl;
  syst_analysis_sep(v_ave_an_Deltamud_col);
  for(size_t i=0;i<8;i++)
    cout<<"an_Deltamud_col: "<<v_ave_an_Deltamud_col[i]<<endl;

  v_ave_an_Deltamud_col_ind=ave_analyses(Deltamud_col_ind);
  cout<<"Deltamud_col_ind: "<<stat_analysis(v_ave_an_Deltamud_col_ind)<<" "<<syst_analysis(v_ave_an_Deltamud_col_ind)<<endl;
  syst_analysis_sep(v_ave_an_Deltamud_col_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_Deltamud_col_ind: "<<v_ave_an_Deltamud_col_ind[i]<<endl;

  v_ave_an_ratio_mu_md_col=ave_analyses(ratio_mu_md_col);
  cout<<"mu/md_col: "<<stat_analysis(v_ave_an_ratio_mu_md_col)<<" "<<syst_analysis(v_ave_an_ratio_mu_md_col)<<endl;
  syst_analysis_sep(v_ave_an_ratio_mu_md_col);
  for(size_t i=0;i<8;i++)
    cout<<"an_ratio_mu_md_col: "<<v_ave_an_ratio_mu_md_col[i]<<endl;

  v_ave_an_ratio_mu_md_col_ind=ave_analyses(ratio_mu_md_col_ind);
  cout<<"mu/md_col_ind: "<<stat_analysis(v_ave_an_ratio_mu_md_col_ind)<<" "<<syst_analysis(v_ave_an_ratio_mu_md_col_ind)<<endl;
  syst_analysis_sep(v_ave_an_ratio_mu_md_col_ind);
  for(size_t i=0;i<8;i++)
    cout<<"an_ratio_mu_md_col_ind: "<<v_ave_an_ratio_mu_md_col_ind[i]<<endl;

  v_ave_an_mu=ave_analyses(mu_q);
  cout<<"mu: "<<stat_analysis(v_ave_an_mu)<<" "<<syst_analysis(v_ave_an_mu)<<endl;
  syst_analysis_sep(v_ave_an_mu);
  for(size_t i=0;i<8;i++)
    cout<<"an_mu: "<<v_ave_an_mu[i]<<endl;

  v_ave_an_md=ave_analyses(md_q);
  cout<<"md: "<<stat_analysis(v_ave_an_md)<<" "<<syst_analysis(v_ave_an_md)<<endl;
  syst_analysis_sep(v_ave_an_md);
  for(size_t i=0;i<8;i++)
    cout<<"an_md: "<<v_ave_an_md[i]<<endl;

  ////////////////////D meson////////////////////
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t daM2D_QCD(raw_data.size());
	dbvec_t daM2D_QCD_ind(raw_data.size());
      
	for(size_t iens=0;iens<raw_data.size();iens++)
         {
           daM2D_QCD[iens]=4.0*Deltamud[ind_an({input_an_id,an_flag})]*MD_sl_s[iens+input_an_id*raw_data.size()]*aMD[iens+input_an_id*raw_data.size()]/1000.0;
           daM2D_QCD_ind[iens]=4.0*Deltamud_ind[ind_an({input_an_id,an_flag})]*MD_sl_s[iens+input_an_id*raw_data.size()]*aMD[iens+input_an_id*raw_data.size()]/1000.0;
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

	cout<<"                         QCD D "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dM2D_QCD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_dM2D_QCD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							 daM2D_QCD[iens],daM2D_QCD[iens]));
	
	output_dM2D_QCD[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2D_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2D_QCD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QCD_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);
	
	cout<<"-----------------------------------------------"<<endl;
	
	cout<<"                         QCD D_ind "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_dM2D_QCD_ind;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL1)
	    data_dM2D_QCD_ind.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
							     daM2D_QCD_ind[iens],daM2D_QCD_ind[iens]));
	output_dM2D_QCD_ind[ind_an({input_an_id,an_flag})]=cont_chir_fit_dM2D_QCD(alist,zlist,lat_par[input_an_id].f0,lat_par[input_an_id].B0,data_dM2D_QCD_ind,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_dM2D_QCD_ind_flag%zu_an%zu.xmg",an_flag,input_an_id),an_flag,cov_flag,beta_list);

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
       dM2D[ind_an({input_an_id,an_flag})]=output_dM2D_QCD[ind_an({input_an_id,an_flag})]+output_dM2D_QED[ind_an({input_an_id,an_flag})];
       dM2D_ind[ind_an({input_an_id,an_flag})]=output_dM2D_QCD_ind[ind_an({input_an_id,an_flag})]+output_dM2D_QED[ind_an({input_an_id,an_flag})];
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
	dbvec_t ml(raw_data.size());
	dbvec_t mc(raw_data.size());
	dbvec_t saMD(raw_data.size());
	dbvec_t FVE_saMD(raw_data.size());
	
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    size_t ibeta=raw_data[iens].ibeta;
	    
	    dboot_t a=1.0/lat_par[input_an_id].ainv[ibeta];
	    ml[iens]=raw_data[iens].aml/lat_par[input_an_id].Z[ibeta]/a;
	    mc[iens]=raw_data[iens].amc/lat_par[input_an_id].Z[ibeta]/a;
	    
	    saMD[iens]=-(eu+ed)*eu*e2*aMD_sl_exch[iens+input_an_id*raw_data.size()]-2.0*sqr(eu)*e2*aMD_sl_selftad_revins[iens+input_an_id*raw_data.size()]-(sqr(eu)+sqr(ed))*e2*aMD_sl_selftad[iens+input_an_id*raw_data.size()]-2.0*aDeltam_cr_u[iens+input_an_id*raw_data.size()]*MD_sl_p_revins[iens+input_an_id*raw_data.size()]-(aDeltam_cr_u[iens+input_an_id*raw_data.size()]+aDeltam_cr_d[iens+input_an_id*raw_data.size()])*MD_sl_p[iens+input_an_id*raw_data.size()]+2.0*a*mc[iens]*MD_sl_s_revins[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(eu))+2.0*a*ml[iens]*MD_sl_s[iens+input_an_id*raw_data.size()]/(Z_QED[iens+input_an_id*raw_data.size()]*2.0/(sqr(ed)+sqr(eu)));
	    FVE_saMD[iens]=FVE_M(aMD[iens+input_an_id*raw_data.size()],raw_data[iens].L);
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

	cout<<"                         SMD "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_sMD;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL2)
	    data_sMD.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMD[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						    saMD[iens]-FVE_saMD[iens],saMD[iens]));
      
	output_sMD[ind_an({input_an_id,an_flag})]=cont_chir_linear_fit_dM(alist,zlist,data_sMD,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_sMD_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),"$$M_{D^+}+M_{D^0} [GeV]",1.0,0.0,an_flag,1,cov_flag,beta_list);

      
	cout<<"-----------------------------------------------"<<endl;
      }

  v_ave_an_sMD=ave_analyses(output_sMD);
  cout<<"M_{D^+}+M_{D^0}: "<<stat_analysis(v_ave_an_sMD)<<" "<<syst_analysis(v_ave_an_sMD)<<endl;
  syst_analysis_sep(v_ave_an_sMD);
  for(size_t i=0;i<8;i++)
    cout<<"an_sMD: "<<v_ave_an_sMD[i]<<endl;

  //////////////////////////Ds meson////////////////////////
  for(size_t an_flag=0;an_flag<nan_syst;an_flag++)
    for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
      {
	dbvec_t ms(raw_data.size());
	dbvec_t mc(raw_data.size());
	dbvec_t aMDsg(raw_data.size());
	dbvec_t FVE_aMDsg(raw_data.size());
	
	for(size_t iens=0;iens<raw_data.size();iens++)
	  {
	    size_t ens_id=raw_data[iens].iult;
	    size_t ibeta=raw_data[iens].ibeta;
	    
	    dboot_t a=1.0/lat_par[input_an_id].ainv[ibeta];
	    
	    boot_init_t &bi=jack_index[input_an_id][ens_id];

	    ms[iens]=raw_data[iens].ams/lat_par[input_an_id].Z[ibeta]/a;
	    mc[iens]=raw_data[iens].amc/lat_par[input_an_id].Z[ibeta]/a;

	    dboot_t aMDs_sl_exch=dboot_t(bi,raw_data[iens].Ds_SL_exch);
	    dboot_t aMDs_sl_selftad=dboot_t(bi,raw_data[iens].Ds_SL_selftad);
	    dboot_t MDs_sl_p=dboot_t(bi,raw_data[iens].Ds_SL_p);
	    dboot_t aMDs_sl_selftad_revins=dboot_t(bi,raw_data[iens].Ds_SL_selftad_revins);
	    dboot_t MDs_sl_p_revins=dboot_t(bi,raw_data[iens].Ds_SL_p_revins);
	    dboot_t MDs_sl_s=dboot_t(bi,raw_data[iens].Ds_SL_s);
	    dboot_t MDs_sl_s_revins=dboot_t(bi,raw_data[iens].Ds_SL_s_revins);
	    
	    aMDsg[iens]=-ed*eu*e2*aMDs_sl_exch-sqr(eu)*e2*aMDs_sl_selftad_revins-aDeltam_cr_u[iens+input_an_id*raw_data.size()]*MDs_sl_p_revins-aDeltam_cr_d[iens+input_an_id*raw_data.size()]*MDs_sl_p-sqr(ed)*e2*aMDs_sl_selftad+a*mc[iens]*MDs_sl_s_revins/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(eu))+ms[iens]*a*MDs_sl_s/(Z_QED[iens+input_an_id*raw_data.size()]/sqr(ed));
	    FVE_aMDsg[iens]=FVE_M(aMDs[iens+input_an_id*raw_data.size()],raw_data[iens].L);
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

	cout<<"                         MDs "<<endl;
	cout<<endl;

	vector<cont_chir_fit_data_t> data_MDs;
	for(size_t iens=0;iens<raw_data.size();iens++)
	  if(raw_data[iens].useforL2)
	    data_MDs.push_back(cont_chir_fit_data_t(raw_data[iens].aml,raw_data[iens].ams,aMDs[iens+input_an_id*raw_data.size()],raw_data[iens].ibeta,raw_data[iens].L,
						    aMDsg[iens]-FVE_aMDsg[iens],aMDsg[iens]));
      
	output_MDs[ind_an({input_an_id,an_flag})]=cont_chir_linear_fit_dM(alist,zlist,data_MDs,lat_par[input_an_id].ml,combine("plots/cont_chir_fit_MDs_flag%zu_an%zu_sub%s.xmg",an_flag,input_an_id,"%s"),"$$M_{Ds^+} [GeV]",1.0,0.0,an_flag,1,cov_flag,beta_list);

      
	cout<<"-----------------------------------------------"<<endl;
      }

  v_ave_an_MDs=ave_analyses(output_MDs);
  cout<<"M_{Ds^+}: "<<stat_analysis(v_ave_an_MDs)<<" "<<syst_analysis(v_ave_an_MDs)<<endl;
  syst_analysis_sep(v_ave_an_MDs);
  for(size_t i=0;i<8;i++)
    cout<<"an_MDs: "<<v_ave_an_MDs[i]<<endl;
  */
  return 0;
}
