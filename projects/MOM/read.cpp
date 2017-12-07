#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#define EXTERN_READ
 #include <read.hpp>

#include <corrections.hpp>
#include <evolutions.hpp>
#include <ingredients.hpp>
#include <prop.hpp>

void read_input(const string &input_path)
{
  using namespace glb;
  
  raw_file_t input(input_path,"r");
  
  const size_t Ls=input.read<size_t>("L"); //!< lattice spatial size
  L[0]=input.read<size_t>("T");
  
  const string act_str=input.read<string>("Action"); //!< action name
  auto act_key=gaz::decr.find(act_str); //!< key in the map of act
  if(act_key==gaz::decr.end()) CRASH("Unable to decript %s",act_str.c_str());
  act=act_key->second;
  
  beta=input.read<double>("Beta");
  plaq=input.read<double>("Plaq");
  nm=input.read<double>("Nm");
  am.resize(nm);
  for(size_t im=0;im<nm;im++) am[im]=input.read<double>();
  am_max=*max_element(am.begin(),am.end())*1.1;
  am_min=*min_element(am.begin(),am.end())/1.1;
  nr=input.read<double>("Nr");
  
  im_sea=input.read<double>("ImSea");
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  set_njacks(input.read<size_t>("NJacks")); //!< number of jacks
  
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  nhits=input.read<size_t>("NHits");
  nhits_to_use=input.read<size_t>("NHitsToUse");
  
  tmin=input.read<size_t>("Tmin");
  tmax=input.read<size_t>("Tmax");
  
  use_QED=input.read<bool>("UseQED");
  print_each_mom=input.read<int>("PrintEachMom");
  
  ainv=input.read<double>("aInv");
  Nf=ev::Nf_t_of_Nf(input.read<int>("Nf"));
  
  ord=input.read<size_t>("Ord");
  
  //////////////////////////////////////////////////
  
  //! sufffix if a number of hits is different from 1
  if(nhits>1) suff_hit="_hit_%zu";
  
  g2=6.0/beta;
  g2tilde=g2/plaq;
  cout<<"g2tilde: "<<g2tilde<<endl;
  
  //set spatial sizes
  V=L[0]*pow(Ls,NDIM-1);
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  //initialize momenta
  set_glb_list_of_moms(mom_list_path);
}

void prepare_list_of_confs()
{
  using namespace glb;
  
  //create the list of all confs available
  string test_path="out/%04zu/fft_S_M0_R0_0";
  if(nhits>1) test_path+="_hit_0";
  conf_list=get_existing_paths_in_range(test_path,conf_range);
  if(conf_list.size()==0) CRASH("list of configurations is empty! check %s ",test_path.c_str());
  clust_size=trim_to_njacks_multiple(conf_list,true);
  
  conf_ind.set_ranges({{"ijack",njacks},{"i_in_clust",clust_size}});
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  i_in_clust_ihit_ind.set_ranges({{"i_in_clust",clust_size},{"ihit",nhits_to_use}});
  im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  r_r_iZbil_ind.set_ranges({{"r",nr},{"r",nr},{"iZbil",nZbil}});
  im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  im_r_ijackp1_ind=im_r_ind*index_t({{"ijack",njacks+1}});
}

namespace
{
  using namespace glb;
  index_t get_im_r_iconf_ihit_ikind_ind()
  {
    return im_r_ind*index_t({{"conf",conf_list.size()},{"hit",nhits},{"kind",m_r_mom_conf_props_t::nprop_kind()}});
  }
}

vector<raw_file_t> setup_read_all_props_mom(const vector<size_t> &conf_list)
{
  using namespace glb;
  
  const index_t im_r_iconf_ihit_ikind_ind=get_im_r_iconf_ihit_ikind_ind();
  vector<raw_file_t> files(im_r_iconf_ihit_ikind_ind.max());
  
#pragma omp parallel for
  for(size_t i=0;i<im_r_iconf_ihit_ikind_ind.max();i++)
    {
      const vector<size_t> comps=im_r_iconf_ihit_ikind_ind(i);
      const size_t im=comps[0];
      const size_t r=comps[1];
      const size_t iconf=comps[2];
      const size_t ihit=comps[3];
      const size_t ikind=comps[4];
      const string path_base=combine("out/%04zu/fft_",conf_list[iconf]);
      const string path_suff=combine(suff_hit.c_str(),ihit);
      
      files[i].open(path_base+get_prop_tag(im,r,ikind)+path_suff,"r");
    }
  
  return files;
}

vector<m_r_mom_conf_props_t> read_all_props_mom(vector<raw_file_t> &files,const size_t i_in_clust_ihit,const size_t imom)
{
  //! output
  vector<m_r_mom_conf_props_t> props(im_r_ijack_ind.max());
  
  //! getting the correct file
  const index_t im_r_iconf_ihit_ikind_ind=get_im_r_iconf_ihit_ikind_ind();
  
  //decompose the outer index
  const vector<size_t> i_in_clust_ihit_comp=i_in_clust_ihit_ind(i_in_clust_ihit);
  const size_t i_in_clust=i_in_clust_ihit_comp[0],ihit=i_in_clust_ihit_comp[1];
  
  //! index of all that must be read
  const index_t im_r_ijack_ikind_ind=im_r_ijack_ind*index_t({{"ikind",m_r_mom_conf_props_t::nprop_kind()}});
#pragma omp parallel for
  for(size_t im_r_ijack_ikind=0;im_r_ijack_ikind<im_r_ijack_ikind_ind.max();im_r_ijack_ikind++)
    {
      const vector<size_t> im_r_ijack_ikind_comps=im_r_ijack_ikind_ind(im_r_ijack_ikind);
      const size_t im=im_r_ijack_ikind_comps[0];
      const size_t r=im_r_ijack_ikind_comps[1];
      const size_t ijack=im_r_ijack_ikind_comps[2];
      const size_t ikind=im_r_ijack_ikind_comps[3];
      
      //! index of the conf built from ijack and i_in_clust
      const size_t iconf=i_in_clust+clust_size*ijack;
      
      //! index of the file to use
      const size_t im_r_iconf_ihit_ikind=im_r_iconf_ihit_ikind_ind({im,r,iconf,ihit,ikind});
      
      //! index of the output propagator
      const size_t im_r_ijack=im_r_ijack_ind({im,r,ijack});
      
      read_prop(props[im_r_ijack].kind[ikind],files[im_r_iconf_ihit_ikind],m_r_mom_conf_props_t::coeff_to_read(ikind,r),imom);
    }
  
  return props;
}
