#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>
#include <unistd.h>

// //class to minimze
// template <class T>
// class plan_minimizer_t : public minimizer_fun_t
// {
//   const plan_fit_data_t<T>& data;
//   const size_t &iel;
  
// public:
  
//   plan_minimizer_t(const plan_fit_data_t<T>& data,const size_t& iel) : data(data),iel(iel) {}
  
//   //! compute the function
//   double operator()(const vector<double> &p) const
//   {
//     return plan_chi2(data,p,iel);
//   }
  
//   double Up() const {return 1;}
// };

int main(int narg,char **arg)
{
  set_njacks(100);
  
  if(narg<3) CRASH("Use: %s file nx",arg[0]);
  string path=arg[1];
  size_t nx=(size_t)atoi(arg[2]);
  
  //open the file
  ifstream file(path);
  if(not file.good()) CRASH("Unable to open %s",arg[1]);
  
  //read
  plan_fit_data_t<djack_t> data;
  int seed=0;
  while(not file.eof())
    {
      vector<double> x(nx+1);
      x[0]=1.0;
      double y,ey;
      for(size_t i=1;i<=nx;i++) file>>x[i];
      file>>y>>ey;
      
      djack_t j;
      j.fill_gauss(y,ey,seed++);
      if(not file.eof())
	data.push_back(make_tuple(x,j));
    }
  
  //compute data size
  const size_t npoints=data.size();
  cout<<"NPoints: "<<npoints<<endl;
  
  const djvec_t res=plan_fit(data);
  
  //const double chi2=plan_chi2(data,res,0);
  //cout<<"Chi2: "<<chi2<<" / "<<npoints-(nx+1)<<endl;
  
  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
  for(size_t i=0;i<=nx;i++)
    cout<<res[i]<<endl;
  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
  
  minimizer_pars_t pars;
  for(size_t i=0;i<=nx;i++)
    pars.add(combine("c%zu",i),0.0,1.0);
  
  // size_t iel;
  // c_t<djack_t> c(data,iel);
  // minimizer_t miner(c,pars);
  
  // djvec_t res_miner(nx+1);
  // for(iel=0;iel<=njacks;iel++)
  //   {
  //     vector<double> res_miner_per_jack=miner.minimize();
  //     for(size_t i=0;i<=nx;i++)
  // 	res_miner[i][iel]=res_miner_per_jack[i];
  //   }
  
  // cout<<res_miner.ave_err()<<endl;
  
  return 0;
}
