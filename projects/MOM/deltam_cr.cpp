#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_DELTAM_CR
 #include <deltam_cr.hpp>

#include <contractions.hpp>
#include <prop.hpp>
#include <read.hpp>
#include <timings.hpp>

void get_deltam_cr()
{
  using namespace glb;
  
  deltam_cr.resize(nm);
  
  //if file exists open it, otherwise compute it
  const string deltam_cr_path="deltam_cr.dat";
  if(file_exists(deltam_cr_path)) deltam_cr.bin_read(deltam_cr_path);
  else
    {
      prepare_list_of_confs();
      
      deltam_cr_time.start();
      for(size_t im=0;im<nm;im++)
	{
	  deltam_cr[im]=compute_deltam_cr(conf_list,tmin,tmax,im,nr,use_QED);
	  cout<<"Deltam cr["<<im<<"]: "<<deltam_cr[im]<<endl;
	}
      deltam_cr_time.stop();
      deltam_cr.bin_write(deltam_cr_path);
    }
}

