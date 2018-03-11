#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/deltam_cr.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

perens_t& perens_t::get_meson_mass()
{
  meson_mass.resize(im_im_ind.max());
  
  //if file exists open it, otherwise compute it
  const string meson_mass_path=dir_path+"/meson_mass.dat";
  if(file_exists(meson_mass_path))
    {
      cout<<"File "<<meson_mass_path<<" found, opening"<<endl;
      deltam_cr.bin_read(meson_mass_path);
    }
  else
    {
      cout<<"File "<<meson_mass_path<<" not found, computing"<<endl;
      prepare_list_of_confs();
      
      meson_mass_time.start();
      for(size_t im1=0;im1<nm;im1++)
	for(size_t im2=0;im2<nm;im2++)
	  meson_mass[im_im_ind({im1,im2})]=compute_meson_mass(im1,im2);
      meson_mass_time.stop();
      meson_mass.bin_write(meson_mass_path);
    }
  
  grace_file_t meson_mass2_plot(dir_path+"/plots/M2.xmg");
  meson_mass2_plot.new_data_set();
  for(size_t im1=0;im1<nm;im1++)
    for(size_t im2=im1;im2<nm;im2++)
      meson_mass2_plot.write_ave_err(am[im1]+am[im2],sqr(meson_mass[im_im_ind({im1,im2})]).ave_err());
  
  return *this;
}
