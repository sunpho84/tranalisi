#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

size_t noa=8,nbeta=3,nens=15,nboots=100;

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
  
  djack_t pi_mass,pi_SL;
  //k_mass,k_SL,D_mass,D_SL;
  //djack_t deltam_cr;
};

vector<lat_par_t> lat_par(noa);
vector<ens_data_t> raw_data;

int main(int narg,char **arg)
{
  set_njacks(15);
  
  //open input file
  string name="input_global_test.txt";
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
      
      //store in the raw_data vector
      raw_data.push_back(temp);
    }
  
  //test the inputx
  cout<<raw_data[0].pi_mass.ave_err()<<endl;

  raw_file_t file(ens_pars,"r");
  
  boot_init_t jack_index[noa][nens];

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
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      for(size_t iboot=0;iboot<nboots;iboot++)
	file.read(lat_par[ia].ainv[ibeta][iboot]);
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ia=0;ia<noa/2;ia++)
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      for(size_t iboot=0;iboot<nboots;iboot++)
	file.read(lat_par[ia].Z[ibeta][iboot]);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t ia=noa/2;ia<noa;ia++)
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      for(size_t iboot=0;iboot<nboots;iboot++)
	file.read(lat_par[ia].Z[ibeta][iboot]);
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t ia=0;ia<noa;ia++)
    for(size_t iens=0;iens<nens;iens++)
      for(size_t iboot=0;iboot<nboots;iboot++)
	file.read(jack_index[ia][iens][iboot]);
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].dB0=read_boot(file);
  
  // for(size_t iens=0;iens<nens;iens++)
  //   cout<<jack_index[0][iens][0]<<endl;

  //cout<<lat_par[0].ml.ave_err()<<endl;
  
  return 0;
}
