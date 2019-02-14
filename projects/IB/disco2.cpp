#include <tranalisi.hpp>

int T,L;

vector<complex<double>> read_vector(const string &path,int n)
{
  raw_file_t data(path,"r");
  
  vector<complex<double>> d;
  
  for(int i=0;i<n;i++)
    d.emplace_back(data.read<double>(),data.read<double>());
  
  return d;
}

int main(int narg,char **arg)
{
  string name="input.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  T=input.read<int>("T");
  L=input.read<int>("L");
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  set_njacks(input.read<int>("NJacks"));
  
  const int nposs_confs=(conf_range.end-conf_range.start)/conf_range.each;
  cout<<"Npossible confs: "<<nposs_confs<<endl;
  const int clust_size=nposs_confs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  const int nconfs=clust_size*njacks;
  cout<<"Nconfs: "<<nconfs<<endl;
  
  const int nhits=input.read<int>("NHits");
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      const int conf=iconf*conf_range.each+conf_range.start;
      vector<complex<double>> EU1_stoch=read_vector(combine("%04d/EU1_stoch",conf),nhits);
      vector<complex<double>> EU2_stoch=read_vector(combine("%04d/EU2_stoch",conf),nhits);
      vector<complex<double>> EU4_stoch=read_vector(combine("%04d/EU4_stoch",conf),nhits);
      vector<complex<double>> EU5_stoch=read_vector(combine("%04d/EU5_stoch",conf),nhits*(nhits-1)/2);
      vector<complex<double>> EU6_stoch=read_vector(combine("%04d/EU6_stoch",conf),nhits*(nhits-1)/2);
      
    }
  
  return 0;
}
