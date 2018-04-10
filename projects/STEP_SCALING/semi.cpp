#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=32;
double m_list[]={0.10,0.20,0.30};
double th_list[]={0.00,-0.10,0.10,-0.20,0.20};
const size_t nth=5,nm=3,nr=2;
const index_t i2pts_ind({{"m1",nm},{"r1",nr},{"m0",nm},{"r0",nr},{"th1",nth},{"th0",nth},{"ri",2},{"T",T}});
const index_t i3pts_ind({{"ms",nm},{"rs",nr},{"m1",nm},{"r1",nr},{"m0",nm},{"r0",nr},{"th1",nth},{"th0",nth},{"ri",2},{"T",T}});
const string base_run="/marconi_work/INF17_lqcd123_0/sanfo/STEP_SCALING/L16_T32_beta6.61/Semileptonic/";

djvec_t retrive_data(const string &outfile,const string &template_path,const index_t &ind)
{
  const range_t range{900,100,40800};
  const size_t ntot_col=2;
  const vector<size_t> cols={0,1};
  
  cout<<"npts: "<<ind.max()<<endl;
  
  djvec_t data;
  if(file_exists(outfile))
    {
      cout<<"File exists, loading"<<endl;
      
      data.resize(ind.max());
      data.bin_read(outfile);
    }
  else
    {
      cout<<"File does not exist, reading confs"<<endl;
      
      data=read_conf_set_t(template_path,range,ntot_col,cols);
      if(data.size()!=ind.max()) CRASH("loaded txt data has size %zu, expected %zu",data.size(),ind.max());
      data.bin_write(outfile);
    }
  
  return data;
}

int main()
{
  set_njacks(15);
  
  cout<<"n3pts: "<<i3pts_ind.max()<<endl;
  
  
  djvec_t data_2pts=retrive_data("2pts.dat",base_run+"out/%05d/mes_contr_2pts",i2pts_ind);
  djvec_t data_3pts=retrive_data("3pts.dat",base_run+"out/%05d/mes_contr_3pts",i3pts_ind);
  
  return 0;
}
