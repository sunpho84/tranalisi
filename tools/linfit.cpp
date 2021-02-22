#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  set_njacks(1000);
  
  if(narg<2) CRASH("Use: %s file [plot] [min|-] [max|-]",arg[0]);
  
  vector<tuple<double,double,double>> tmp;
  
  ifstream file(arg[1]);
  if(not file.good()) CRASH("Unable to open %s",arg[1]);
  
  vector<double> x;
  vec_ave_err_t _y;
  
  double xi,yi,ei;
  while(file>>xi>>yi>>ei)
    {
      x.push_back(xi);
      _y.push_back({yi,ei});
    }
  
  djvec_t y(_y.size());
  for(size_t i=0;i<_y.size();i++)
    y[i].fill_gauss(_y[i],i+123124);
  
  const auto minMax=minmax_element(x.begin(),x.end());
  const string plotPath=(narg<3)?"":arg[2];
  
  const double min=(narg<4 or (string)arg[3]=="-")?*minMax.first*0.9:strtod(arg[3],NULL);
  const double max=(narg<5 or (string)arg[4]=="-")?*minMax.second*1.1:strtod(arg[4],NULL);
  const djvec_t par=poly_fit(x,y,1,min,max,plotPath);
  
  cout<<par.ave_err()<<endl;
  
  return 0;
}
