#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int main()
{
  set_njacks(15);
  
  const range_t range{900,100,40800};
  const size_t ntot_col=2;
  const vector<size_t> cols={0,1};
  const string template_path_3pts="out/%05d/mes_contr_3pts";
  read_conf_set_t(template_path_3pts,range,ntot_col,cols);
  
  return 0;
}
