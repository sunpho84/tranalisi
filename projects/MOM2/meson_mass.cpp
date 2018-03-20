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
      meson_mass.bin_read(meson_mass_path);
    }
  else
    {
      cout<<"File "<<meson_mass_path<<" not found, computing"<<endl;
      
      prepare_list_of_confs();
      meson_mass_time.start();
      for(size_t im1=0;im1<nm;im1++)
	for(size_t im2=0;im2<nm;im2++)
	  meson_mass[im_im_ind({im1,im2})]=compute_meson_mass(to_string(im1),to_string(im2));
      meson_mass_time.stop();
      meson_mass.bin_write(meson_mass_path);
    }
  
  grace_file_t meson_mass2_plot(dir_path+"/plots/M2.xmg");
  meson_mass2_plot.new_data_set();
  for(size_t im1=0;im1<nm;im1++)
    for(size_t im2=im1;im2<nm;im2++)
      meson_mass2_plot.write_ave_err(am[im1]+am[im2],sqr(meson_mass[im_im_ind({im1,im2})]).ave_err());
  
  //sea meson
  const string meson_mass_sea_path=dir_path+"/meson_mass_sea.dat";
  if(im_sea>=0 and im_sea<(int)nm) meson_mass_sea=meson_mass[im_im_ind({(size_t)im_sea,(size_t)im_sea})];
  else
    if(file_exists(meson_mass_sea_path))
      {
	cout<<"File "<<meson_mass_sea_path<<" found, opening"<<endl;
	meson_mass_sea.bin_read(meson_mass_sea_path);
      }
    else
      {
	cout<<"File "<<meson_mass_sea_path<<" not found, computing"<<endl;
	
	prepare_list_of_confs();
	meson_mass_time.start();
	meson_mass_sea=compute_meson_mass("sea","sea");
	meson_mass_time.stop();
	meson_mass_sea.bin_write(meson_mass_sea_path);
      }
  cout<<"Sea meson mass: "<<smart_print(meson_mass_sea.ave_err())<<endl;
  
  return *this;
}
