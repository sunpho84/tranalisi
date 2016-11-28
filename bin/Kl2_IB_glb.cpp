#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

class ens_data_t
{
public:
  size_t iult,ibeta,L;
  double aml,ams,amc;
  string path;
  
  djack_t pi_mass,pi_SL;
};

vector<ens_data_t> raw_data;

int main(int narg,char **arg)
{
  set_njacks(15);
  
  //open input file
  string name="analysis_input";
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
      raw_file_t obs_pi_file(combine("%s/obs_pi",temp.path.c_str()),"r");
      obs_pi_file.bin_read(temp.pi_mass);
      obs_pi_file.bin_read(temp.pi_SL);
      
      //store in the raw_data vector
      raw_data.push_back(temp);
    }
  
  //test the inputx
  cout<<raw_data[0].pi_mass.ave_err()<<endl;
  
  return 0;
}
