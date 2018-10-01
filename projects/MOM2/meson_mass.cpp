#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <MOM2/deltam.hpp>
#include <MOM2/perens.hpp>
#include <MOM2/timings.hpp>

void perens_t::bin_read_meson_mass()
{
  meson_mass.bin_read(meson_mass_path());
}

void perens_t::bin_write_meson_mass()
{
  meson_mass.bin_write(meson_mass_path());
}

void perens_t::bin_read_meson_mass_QED()
{
  meson_mass_QED.bin_read(meson_mass_QED_path());
}

void perens_t::bin_write_meson_mass_QED()
{
  meson_mass_QED.bin_write(meson_mass_QED_path());
}

void perens_t::bin_read_meson_mass_sea()
{
  meson_mass_sea.bin_read(meson_mass_sea_path());
}

void perens_t::bin_write_meson_mass_sea()
{
  meson_mass_sea.bin_write(meson_mass_sea_path());
}

perens_t& perens_t::get_meson_mass()
{
  //if file exists open it, otherwise compute it
  if(file_exists(meson_mass_path()))
    {
      cout<<"File "<<meson_mass_path()<<" found, opening"<<endl;
      bin_read_meson_mass();
    }
  else
    {
      cout<<"File "<<meson_mass_path()<<" not found, computing"<<endl;
      
      prepare_list_of_confs();
      meson_mass_time.start();
      for(size_t im1=0;im1<nm;im1++)
	for(size_t im2=0;im2<nm;im2++)
	  meson_mass[im_im_ind({im1,im2})]=compute_meson_mass(to_string(im1),to_string(im2));
      meson_mass_time.stop();
      bin_write_meson_mass();
    }
  
  grace_file_t meson_mass2_plot(dir_path+"/plots/M2.xmg");
  meson_mass2_plot.new_data_set();
  for(size_t im1=0;im1<nm;im1++)
    for(size_t im2=im1;im2<nm;im2++)
      meson_mass2_plot.write_ave_err(am[im1]+am[im2],sqr(meson_mass[im_im_ind({im1,im2})]).ave_err());
  
  //if file exists open it, otherwise compute it
  if(pars::use_QED)
    {
      if(file_exists(meson_mass_QED_path()))
	{
	  cout<<"File "<<meson_mass_QED_path()<<" found, opening"<<endl;
	  bin_read_meson_mass_QED();
	}
      else
	{
	  cout<<"File "<<meson_mass_QED_path()<<" not found, computing"<<endl;
	  
	  prepare_list_of_confs();
	  meson_mass_time.start();
	  for(size_t im1=0;im1<nm;im1++)
	    for(size_t im2=0;im2<nm;im2++)
	      meson_mass_QED[im_im_ind({im1,im2})]=compute_meson_mass_QED(im1,im2);
	  meson_mass_time.stop();
	  bin_write_meson_mass_QED();
	}
      
      grace_file_t meson_mass_QED2_plot(dir_path+"/plots/M2_QED.xmg");
      meson_mass_QED2_plot.new_data_set();
      for(size_t im1=0;im1<nm;im1++)
	for(size_t im2=im1;im2<nm;im2++)
	  meson_mass_QED2_plot.write_ave_err(meson_mass[im_im_ind({im1,im2})].ave(),meson_mass_QED[im_im_ind({im1,im2})].ave_err());
    }
  
  //sea meson
  if(im_sea>=0 and im_sea<(int)nm) meson_mass_sea=meson_mass[im_im_ind({(size_t)im_sea,(size_t)im_sea})];
  else
    if(file_exists(meson_mass_sea_path()))
      {
	cout<<"File "<<meson_mass_sea_path()<<" found, opening"<<endl;
	bin_read_meson_mass_sea();
      }
    else
      {
	cout<<"File "<<meson_mass_sea_path()<<" not found, computing"<<endl;
	
	prepare_list_of_confs();
	meson_mass_time.start();
	meson_mass_sea=compute_meson_mass("sea","sea");
	meson_mass_time.stop();
	bin_write_meson_mass_sea();
      }
  cout<<"Sea meson mass: "<<smart_print(meson_mass_sea.ave_err())<<endl;
  
  return *this;
}

vector<perens_t::task_t> perens_t::get_meson_mass_tasks(const vector<const perens_t*> &ens)
{
  vector<const djvec_t*> in_meson_mass;
  
  for(auto &e : ens)
    in_meson_mass.push_back(&e->meson_mass);
  
  vector<task_t> meson_mass_tasks={{&meson_mass,in_meson_mass,im_im_ind,"meson_mass",QED_task}};
  
  return meson_mass_tasks;
}
