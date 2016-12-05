#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

size_t noa=8,nbeta=3,nboots=100;

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
  
  djack_t pi_mass,pi_SL,k_mass,k_SL_exch,k_SL_selftad,k_SL_s,k_SL_p,D_mass,D_SL_exch,D_SL_selftad,D_SL_s,D_SL_p;
  djack_t deltam_cr;
};
/*
template< class T> T fit_delta_pion_square()
{return e2*sqr(f0)*(2*par[0]/pow(f0,4)-(3+4*par[0]/pow(f0,4))*dB0*aml/a/Z*log(dB0*aml/a/Z/4)/sqr(4*M_PI*f0)+par[1]*dB0*aml/a/Z/sqr(4*M_PI*f0))+par[2]*sqr(a);}
*/
vector<lat_par_t> lat_par(noa);
vector<ens_data_t> raw_data;

template <class T1,class T2> T1 FVE_d2M(const T1 &M,const T2 &Lphys)
{
  const double FVE_k=2.837297;
  double alpha=1.0/137;
  return -FVE_k*alpha/Lphys*(M+2.0/Lphys);
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

const bool chir_an_flag(size_t ia)
{return ia%2==0;}

int main(int narg,char **arg)
{
  set_njacks(15);
  
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
      obs_pi_file.bin_read(temp.pi_SL);

      //read the observable of the kaon
      raw_file_t obs_k_file(combine("%s/k_obs",temp.path.c_str()),"r");
      obs_k_file.bin_read(temp.k_mass);
      obs_k_file.bin_read(temp.k_SL_exch);
      obs_k_file.bin_read(temp.k_SL_selftad);
      obs_k_file.bin_read(temp.k_SL_s);
      obs_k_file.bin_read(temp.k_SL_p);

      //read the observable of the D meson
      raw_file_t obs_D_file(combine("%s/D_obs",temp.path.c_str()),"r");
      obs_D_file.bin_read(temp.D_mass);
      obs_D_file.bin_read(temp.D_SL_exch);
      obs_D_file.bin_read(temp.D_SL_selftad);
      obs_D_file.bin_read(temp.D_SL_s);
      obs_D_file.bin_read(temp.D_SL_p);

      //read deltam_cr (ud)
      raw_file_t deltam_cr_file(combine("%s/ud_fit_deltam_cr",temp.path.c_str()),"r");
      deltam_cr_file.bin_read(temp.deltam_cr);
      
      //store in the raw_data vector
      raw_data.push_back(temp);
    }
  
  //test the inputx
  cout<<raw_data[0].pi_mass.ave_err()<<endl;

  raw_file_t file(ens_pars,"r");
  
  boot_init_t jack_index[noa][nens_used];

  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ian=0;ian<noa;ian++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ian].ainv[ibeta][iboot]);
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ian=0;ian<noa/2;ian++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ian].Z[ibeta][iboot]);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ian=noa/2;ian<noa;ian++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ian].Z[ibeta][iboot]);
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t ian=0;ian<noa;ian++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[ian][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t ian=0;ian<noa;ian++) lat_par[ian].B0=read_boot(file)/2.0;
  
  // for(size_t iens=0;iens<nens_used;iens++)
  //   cout<<jack_index[0][iens][0]<<endl;

  //cout<<lat_par[0].ml.ave_err()<<endl;

  //loop
  cout<<endl;
  for(size_t ian=0;ian<noa;ian++)
    {
      dbvec_t ml(raw_data.size());
      dbvec_t a2M2Pi(raw_data.size());
      dbvec_t da2M2Pi(raw_data.size());
      dbvec_t FVE_da2M2Pi(raw_data.size());
      dbvec_t da2M2K_QED(raw_data.size());
      dbvec_t FVE_da2M2K(raw_data.size());
      dbvec_t epsilon_gamma(raw_data.size());
      dbvec_t epsilon_gamma_minusFVE(raw_data.size());
      
      for(size_t iens=0;iens<raw_data.size();iens++)
	{
	  size_t ens_id=raw_data[iens].iult;
	  size_t ibeta=raw_data[iens].ibeta;
	  
	  dboot_t a=1.0/lat_par[ian].ainv[ibeta];
	  dboot_t Lphys=raw_data[iens].L*a;
	  ml[iens]=raw_data[iens].aml/lat_par[ian].Z[ibeta]/a;
	  
	  boot_init_t &bi=jack_index[ian][ens_id];
	  
	  dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*lat_par[ian].Z[ibeta]*(6*log(2*a)-22.596)/(32*sqr(M_PI)));
	  dboot_t Deltam_cr_u=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(eu)/a;
	  dboot_t Deltam_cr_d=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(ed)/a;
	  
	  dboot_t aMPi=dboot_t(bi,raw_data[iens].pi_mass);
	  dboot_t aMPi_sl=dboot_t(bi,raw_data[iens].pi_SL);
	  
	  a2M2Pi[iens]=aMPi*aMPi;
	  da2M2Pi[iens]=aMPi*sqr(eu-ed)*e2*aMPi_sl;
	  FVE_da2M2Pi[iens]=FVE_d2M(aMPi,raw_data[iens].L);
	  
	  dboot_t aMK=dboot_t(bi,raw_data[iens].k_mass);
	  dboot_t aMK_sl_exch=dboot_t(bi,raw_data[iens].k_SL_exch)/a;
	  dboot_t aMK_sl_selftad=dboot_t(bi,raw_data[iens].k_SL_selftad)/a;
	  dboot_t MK_sl_s=dboot_t(bi,raw_data[iens].k_SL_s);
	  dboot_t MK_sl_p=dboot_t(bi,raw_data[iens].k_SL_p);
	  
	  dboot_t daMK_QED=
	    -2*ml[iens]*MK_sl_s/Z_QED
	    -(Deltam_cr_u-Deltam_cr_d)*MK_sl_p
	    +(sqr(eu)-sqr(ed))*e2*(aMK_sl_exch-aMK_sl_selftad);
	  da2M2K_QED[iens]=daMK_QED*2*aMK;
	  FVE_da2M2K[iens]=FVE_d2M(aMK,raw_data[iens].L);
	}
      
      //prepare the list of a and z
      dbvec_t alist(nbeta),zlist(nbeta);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	{
	  alist[ibeta]=1.0/lat_par[ian].ainv[ibeta];
	  zlist[ibeta]=lat_par[ian].Z[ibeta];
	}
      
      //data to fit
      vector<cont_chir_fit_data_t> data_pi;
      //vector<cont_chir_fit_data_t> data_k;
      for(size_t iens=0;iens<raw_data.size();iens++)
	data_pi.push_back(cont_chir_fit_data_t(raw_data[iens].aml,
					    raw_data[iens].ibeta,
					    raw_data[iens].L,
					    da2M2Pi[iens]-FVE_da2M2Pi[iens]));
      
      cont_chir_fit(alist,zlist,lat_par[ian].f0,lat_par[ian].B0,data_pi,lat_par[ian].ml,combine("plots/cont_chir_fit_dM2Pi_an%zu.xmg",ian),chir_an_flag(ian));
      
      //data to fit
      vector<cont_chir_fit_data_t> data_K;
      //vector<cont_chir_fit_data_t> data_k;
      for(size_t iens=0;iens<raw_data.size();iens++)
	data_K.push_back(cont_chir_fit_data_t(raw_data[iens].aml,
					    raw_data[iens].ibeta,
					    raw_data[iens].L,
					    da2M2K_QED[iens]-FVE_da2M2K[iens]));
      
      cont_chir_fit(alist,zlist,lat_par[ian].f0,lat_par[ian].B0,data_K,lat_par[ian].ml,combine("plots/cont_chir_fit_dM2K_an%zu.xmg",ian),chir_an_flag(ian));
      
    }
  
  return 0;
}
