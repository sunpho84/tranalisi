#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/deltam_cr.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

perens_t& perens_t::get_deltam_cr()
{
  deltam_cr.resize(nm);
  
  //if file exists open it, otherwise compute it
  const string deltam_cr_path=dir_path+"/deltam_cr.dat";
  if(file_exists(deltam_cr_path))
    {
      cout<<"File "<<deltam_cr_path<<" found, opening"<<endl;
      deltam_cr.bin_read(deltam_cr_path);
    }
  else
    {
      cout<<"File "<<deltam_cr_path<<" not found, computing"<<endl;
      prepare_list_of_confs();
      
      deltam_cr_time.start();
      for(size_t im=0;im<nm;im++)
	deltam_cr[im]=compute_deltam_cr(im);
      deltam_cr_time.stop();
      deltam_cr.bin_write(deltam_cr_path);
    }
  
  for(size_t im=0;im<nm;im++)
    cout<<"Deltam cr["<<im<<"]: "<<deltam_cr[im]<<endl;
  
  return *this;
}
