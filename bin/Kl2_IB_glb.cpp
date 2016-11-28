#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

class slopemass_t
{
public:
  size_t iult;
  size_t ibeta;
  size_t L;
  double aml;
  double ams;
  double amc;
  string path;
  
  djack_t pi_mass,pi_SL;
  
};

vector<slopemass_t> raw_data;

int main(int narg,char **arg)
{
  set_njacks(15);
  
  string name="analysis_input";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  string path_pars=input.read<string>("UltimatePath");
  size_t nens_used=input.read<int>("NEnsamble");
  
  for(size_t iens=0;iens<nens_used;iens++)
    {
      slopemass_t temp;
            
      temp.iult=input.read<size_t>("Ens");
      temp.ibeta=input.read<size_t>("beta");
      temp.L=input.read<size_t>("L");
      temp.aml=input.read<double>("aml");
      temp.ams=input.read<double>("ams");
      temp.amc=input.read<double>("amc");
      temp.path=input.read<string>("path");
      
      raw_file_t obs_pi_file(combine("%s/obs_pi",temp.path.c_str()),"r");
      obs_pi_file.bin_read(temp.pi_mass);
      obs_pi_file.bin_read(temp.pi_SL);
      
      raw_data.push_back(temp);
    }
  
  cout<<raw_data[0].pi_mass.ave_err()<<endl;
  
  return 0;
}
