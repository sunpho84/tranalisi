#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <iostream>

using namespace std;

int main(int narg,char **arg)
{
  set_njacks(3);
  
  dbvec_t a(3),b(3);
  dboot_t c;
  a=3;
  b=1;
  c=4;

  dbvec_t d;

  uminus(1);
  
  d=-pow(log(b/(a+c)),2)+2;
  //has_iterator<boot_t<double>>::value;

  cout<<d[0][1]<<endl;

  ave_err_t av=d[0].ave_err();
  cout<<av<<endl;
  
  //cout<<has_method_size<boot_t<double>>::value<<endl;
  
  //decltype(boot_t<double>::iterator(),0);
  
  dboot_t test1(gauss_filler_t(3,2,124142));
  dboot_t test2(gauss_filler_t(3,2,124142));
  cout<<(test1/test2).ave_err()<<endl;
  cout<<test2.ave_err()<<endl;

  boot_init_t binit;binit.fill(23252523);
  vector<vector<double>> ciccio={{3,4,1},{4,5,4}};
  
  dbvec_t bciccio=bvec_from_jvec(binit,djvec_t(ciccio));
  
  cout<<bciccio.ave_err()<<endl;
  
  //obs_file_t file("/tmp/cicc",2);
  //file.set_col_view({0,1});
  //auto i=file.read(4);
  
  //cout<<endl<<i<<endl;
  
  //djvec_t te(i);

  string template_path="cicc%d";
  range_t file_range={0,1,2};
  size_t ntot_cols=2;
  vector<size_t> cols={0,1};
  size_t nlines=3;

  djvec_t jvec=read_conf_set_t(template_path,file_range,ntot_cols,cols,nlines);
  raw_file_t outf("/tmp/test","w");
  jvec.bin_write(outf);
  
  filter_t(2)(jvec).bin_write("/tmp/202");
  
  //cout<<dbvec_t(binit,jvec).ave_err()<<endl;
  
  //auto fun=[](double x){return sin(0.23561634*x+0.36124252);};
  auto fun=[](double x){return (x+3)*sqr(x-1);};
  cout.precision(16);
  cout<<"======================="<<endl;
  cout<<fun(Brent_solve(fun,-1.5))<<endl;

  index_t<4> ind;
  ind.set_ranges({0,9,2,3});
  cout<<ind({0,9,2,3})<<endl;
  
  //int o=file.read<int>("o");
  //cout<<o<<endl;
  return 0;
}
