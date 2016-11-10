#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <file.hpp>
#include <jack.hpp>
#include <oper.hpp>

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
  
  dbvec_t bciccio(binit,ciccio);
  
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
  
  jvec.filter(2).bin_write("/tmp/202");
  
  cout<<dbvec_t(binit,jvec).ave_err()<<endl;
  
  //int o=file.read<int>("o");
  //cout<<o<<endl;
  return 0;
}
