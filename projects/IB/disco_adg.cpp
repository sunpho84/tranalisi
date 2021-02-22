#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

int main(int narg,char **arg)
{
  set_njacks(1000);
  
  if(narg<3) CRASH("Use: %s file ms [plot] [min|-] [max|-]",arg[0]);
  
  double ms=strtod(arg[2],nullptr);
  
  vector<tuple<double,double,double>> tmp;
  
  ifstream file(arg[1]);
  if(not file.good()) CRASH("Unable to open %s",arg[1]);
  
  vector<double> x,xsh;
  vec_ave_err_t _y;
  
  double xi,yi,ei;
  while(file>>xi>>yi>>ei)
    {
      x.push_back(xi);
      xsh.push_back(xi-ms);
      _y.push_back({yi,ei});
    }
  
  djvec_t y(_y.size()),z(_y.size());
  for(size_t i=0;i<_y.size();i++)
    {
      y[i].fill_gauss(_y[i],i+123124);
      z[i]=y[i]/(x[i]-ms);
    }
  
  const auto minMax=minmax_element(x.begin(),x.end());
  const string plotPath=(narg<4)?"":arg[3];
  
  const double min=(narg<5 or (string)arg[4]=="-")?*minMax.first*0.9:strtod(arg[4],NULL);
  const double max=(narg<6 or (string)arg[5]=="-")?*minMax.second*1.1:strtod(arg[5],NULL);
  const djvec_t pars=poly_fit(xsh,z,1);
  djvec_t pars_true(3);
  pars_true[0]=0.0;
  pars_true[1]=pars[0];
  pars_true[2]=pars[1];
  
  if(plotPath!="")
    {
      grace_file_t plot(plotPath);
      
      plot.write_vec_ave_err(x,y.ave_err(),grace::BLACK,grace::SQUARE);
      plot.write_polygon([&pars_true,ms](const double x)->djack_t{return poly_eval(pars_true,x-ms);},min,max,grace::RED);
      
      const djvec_t pars_par=poly_fit(x,y,2);
      
      plot.write_polygon([&pars_par](const double x)->djack_t{return poly_eval(pars_par,x);},min,max,grace::BLUE);
      
    }
  
  cout<<pars.ave_err()<<endl;
  
  return 0;
}
