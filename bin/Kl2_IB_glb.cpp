#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

size_t noa=8,nbeta=3,nboots=100;

const double e2=4*M_PI/137.;
const double eu=2.0/3;
const double ed=-1.0/3;

const vector<size_t> symbol={grace::SQUARE,grace::CIRCLE,grace::DIAMOND};
const vector<size_t> color={grace::GREEN4,grace::RED,grace::BLUE};

class lat_par_t
{
public:
  dboot_t ml,ms,mc,r0,f0,dB0;
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
{return e2*pow(f0,2)*(2*par[0]/pow(f0,4)-(3+4*par[0]/pow(f0,4))*dB0*aml*ainv[ibeta]/Z[ibeta]*log(dB0*aml*ainv[ibeta]/Z[ibeta]/4)/pow(4*pi*f0,2)+par[1]*dB0*aml*ainv[ibeta]/Z[ibeta]/pow(4*pi*f0,2))+par[2]/pow(ainv[ibeta],2);}
*/
vector<lat_par_t> lat_par(noa);
vector<ens_data_t> raw_data;

template <class T> T FVE_d2M(const T &M,const T &Lphys)
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
      file.no_line();
      file.set_symbol(symbol[ibeta]);
      file.color(color[ibeta]);
      file<<"@type xydy"<<endl;
      for(size_t iens=0;iens<raw_data.size();iens++)
	if(raw_data[iens].ibeta==ibeta)
	  file<<x[iens].ave()<<" "<<y[iens].ave_err()<<endl;
      file.new_set();
    }
}
template <class Tx,class Ty> void plot_ens_data(grace_file_t &&file,const Tx &x,const Ty &y)
{plot_ens_data(file,x,y);}

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
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ia=0;ia<noa;ia++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ia].ainv[ibeta][iboot]);
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ia=0;ia<noa/2;ia++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ia].Z[ibeta][iboot]);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ia=noa/2;ia<noa;ia++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[ia].Z[ibeta][iboot]);
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t ia=0;ia<noa;ia++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_used;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[ia][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].dB0=read_boot(file);
  
  // for(size_t iens=0;iens<nens_used;iens++)
  //   cout<<jack_index[0][iens][0]<<endl;

  //cout<<lat_par[0].ml.ave_err()<<endl;

  //loop
  cout<<endl;
  for(size_t ia=0;ia<noa;ia++)
    {
      dbvec_t ml(raw_data.size());
      dbvec_t M2Pi(raw_data.size());
      dbvec_t dM2Pi(raw_data.size());
      dbvec_t FVE_dM2Pi(raw_data.size());
      dbvec_t dM2K_QED(raw_data.size());
      dbvec_t FVE_dM2K(raw_data.size());
      
      for(size_t iens=0;iens<raw_data.size();iens++)
	{
	  size_t ens_id=raw_data[iens].iult;
	  size_t ibeta=raw_data[iens].ibeta;
	  
	  dboot_t a=1.0/lat_par[ia].ainv[ibeta];
	  dboot_t L=raw_data[iens].L*a;
	  ml[iens]=raw_data[iens].aml/lat_par[ia].Z[ibeta]/a;

	  boot_init_t &bi=jack_index[ia][ens_id];
	  
	  dboot_t Z_QED=1.0/((sqr(ed)-sqr(eu))*e2*lat_par[ia].Z[ibeta]*(6*log(2*a)-22.596)/(32*sqr(M_PI)));
	  dboot_t Deltam_cr_u=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(eu)/a;
	  dboot_t Deltam_cr_d=dboot_t(bi,raw_data[iens].deltam_cr)*e2*sqr(ed)/a;

	  dboot_t MPi=dboot_t(bi,raw_data[iens].pi_mass)/a;
	  dboot_t MPi_sl=dboot_t(bi,raw_data[iens].pi_SL)/a;

	  M2Pi[iens]=MPi*MPi;
	  dM2Pi[iens]=MPi*sqr(eu-ed)*e2*MPi_sl;
	  FVE_dM2Pi[iens]=FVE_d2M(MPi,L);
	  
	  dboot_t MK=dboot_t(bi,raw_data[iens].k_mass)/a;
	  dboot_t MK_sl_exch=dboot_t(bi,raw_data[iens].k_SL_exch)/a;
	  dboot_t MK_sl_selftad=dboot_t(bi,raw_data[iens].k_SL_selftad)/a;
	  dboot_t MK_sl_s=dboot_t(bi,raw_data[iens].k_SL_s);
	  dboot_t MK_sl_p=dboot_t(bi,raw_data[iens].k_SL_p);

	  dboot_t dMK_QED=
	    -2*ml[iens]*MK_sl_s/Z_QED
	    -(Deltam_cr_u-Deltam_cr_d)*MK_sl_p
	    +(sqr(eu)-sqr(ed))*e2*(MK_sl_exch-MK_sl_selftad);
	  dM2K_QED[iens]=dMK_QED*2*MK;
	  FVE_dM2K[iens]=FVE_d2M(MK,L);	  
	}
      
      plot_ens_data(combine("plots/dM2Pi_an%zu.xmg",ia),ml,dM2Pi);
      plot_ens_data(combine("plots/dM2Pi_FVEcorr_an%zu.xmg",ia),ml,dM2Pi-FVE_dM2Pi);
      plot_ens_data(combine("plots/dM2K_QED_an%zu.xmg",ia),ml,dM2K_QED);
      plot_ens_data(combine("plots/dM2K_QED_FVEC_corr_an%zu.xmg",ia),ml,dM2K_QED-FVE_dM2K);
    }

  return 0;
}
