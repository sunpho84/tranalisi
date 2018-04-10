#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const range_t range{900,100,40800};
const size_t ntot_col=2;
const vector<size_t> cols={0,1};

const size_t nth=5,nm=3,nr=2;
const index_t i2pts_ind({{"m1",nm},{"r1",nr},{"m0",nm},{"r0",nr},{"th1",nth},{"th0",nth}});
const index_t i3pts_ind({{"ms",nm},{"rs",nr},{"m1",nm},{"r1",nr},{"m0",nm},{"r0",nr},{"th1",nth},{"th0",nth}});

int main()
{
  set_njacks(15);
  
  cout<<"n2pts: "<<i2pts_ind.max()<<endl;
  cout<<"n3pts: "<<i3pts_ind.max()<<endl;
  
  const string template_path_2pts="out/%05d/mes_contr_2pts";
  djvec_t data_2pts=read_conf_set_t(template_path_2pts,range,ntot_col,cols);
  
  const string template_path_3pts="out/%05d/mes_contr_3pts";
  djvec_t data_3pts=read_conf_set_t(template_path_3pts,range,ntot_col,cols);
  
  
  
  return 0;
}
