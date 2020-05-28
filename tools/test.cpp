#include <tranalisi.hpp>

size_t T,LfrA;
size_t tmin;
size_t ib;
int loc;
range_t file_range;

const array<ave_err_t,3> a_list{ave_err_t{0.0885/0.197,0.0036/0.197},ave_err_t{0.0815/0.197,0.0030/0.197},ave_err_t{0.0619/0.197,0.0018/0.197}};
const array<ave_err_t,3> Zv_list{ave_err_t{0.5920,0.0004},ave_err_t{0.6095,0.0003},ave_err_t{0.6531,0.0002}};

djvec_t read(const string&path,const string& name)
{
  const djvec_t data=read_conf_set_t(path+"/%04d/mes_contr_"+name,file_range,2,{0},T,false);
  const size_t base_nel=T;
  const size_t offset=1;
  const size_t each=2;
  const size_t hw=data.size()/(base_nel*each);
  
  return vec_filter(data,gslice(base_nel*offset,{hw,T},{each*T,1})).symmetrized();
}

djack_t analysis(const string& path)
{
  input_file_t input(path+"/analysis.txt");
  T=input.read<size_t>("T");
  LfrA=input.read<size_t>("L");
  tmin=input.read<size_t>("tmin");
  loc=input.read<size_t>("loc");
  ib=input.read<size_t>("ib");
  input.expect("ConfRange");
  file_range.start=input.read<size_t>();
  file_range.each=input.read<size_t>();
  file_range.end=input.read<size_t>();
  
  string unperturb_path,perturb_path;
  if(loc)
    {
      unperturb_path="M0_R0_0_M0_R0_0";
      perturb_path="M0_R0_F_M0_R0_F";
    }
  else
    {
      unperturb_path="00";
      perturb_path="LL";
    }
  
  djack_t a,Zv;
  
  a.fill_gauss(a_list[ib],ib+235423);
  if(loc)
    Zv.fill_gauss(Zv_list[ib],ib+5686473);
  else
    Zv=1;
  
  const djvec_t LO=read(path,unperturb_path);
  const djvec_t LL=read(path,perturb_path)*sqr(Zv);
  LL.ave_err().write(path+"/LL.xmg");
  const djvec_t R=LL/LO;
  djvec_t dMt=forward_derivative(R);
  
  const djvec_t Mt=effective_mass(LO);
  const djack_t aM=constant_fit(Mt,tmin,T/2,path+"/P5P5.xmg");
  cout<<"M: "<<aM.ave_err()<<endl;
  
  djvec_t C(T/2+1);
  for(size_t t=0;t<T/2;t++)
    C[t]=(T/2-t)*tanh(Mt[t]*(T/2-t));
  
  dMt[0]=dMt[1];
  for(size_t t=1;t<T/2;t++)
    dMt[t]/=C[t+1]-C[t];
  dMt[T/2]=dMt[T/2-1];
  
  const djack_t dM=constant_fit(dMt,tmin,T/2,path+"/dM.xmg");
  cout<<"dM: "<<smart_print(dM)<<endl;
  
  const double e2=4*3.14159/137;
  const double k=2.8373;
  
  const djack_t corr=-dM*aM/sqr(a)*e2;
  cout<<corr.ave_err()<<endl;
  
  const djack_t FVE=e2*k/(4*3.14159)*(aM/LfrA+2/sqr(LfrA))/sqr(a);
  cout<<FVE.ave_err()<<endl;
  
  const djack_t corr_corr=corr+FVE;
  cout<<corr_corr.ave_err()<<endl;
  
  cout<<endl;
  
  return dM;
}

int main(int narg,char **arg)
{
  if(narg<3)
    CRASH("Use %s cons loc",arg[0]);
  
  set_njacks(15);
  
  const djack_t cons=analysis(arg[1]);
  const djack_t loc=analysis(arg[2]);
  
  djack_t a;
  a.fill_gauss(a_list[ib],ib+235423);
  
  const djack_t rat=(cons/loc-1.0)/sqr(a);
  
  cout<<endl;
  cout<<"ratio: "<<endl;
  cout<<smart_print(rat)<<endl;
  
  return 0;
}
