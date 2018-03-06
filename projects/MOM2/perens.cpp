#define EXTERN_PERENS
 #include <MOM2/perens.hpp>

perens_t& perens_t::read_pars(const string &name)
{
  dir_path=name;
  
  raw_file_t input(dir_path+"/pars.txt","r");
  
  const size_t Ls=input.read<size_t>("L"); //!< lattice spatial size
  L[0]=input.read<size_t>("T");
  
  beta=input.read<double>("Beta");
  plaq=input.read<double>("Plaq");
  //load masses and r
  nm=input.read<double>("Nm");
  am.resize(nm);
  for(size_t im=0;im<nm;im++) am[im]=input.read<double>();
  nr=input.read<double>("Nr");
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  set_comp_list_of_moms(dir_path+"/"+mom_list_path,pars::filter_thresh);
  
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  prop_hadr_path=input.read<string>("PropHadrPath");
  prop_lep_path=input.read<string>("PropLepPath");
  
  nhits=input.read<size_t>("NHits");
  nhits_to_use=input.read<size_t>("NHitsToUse");
  
  tmin=input.read<size_t>("Tmin");
  tmax=input.read<size_t>("Tmax");
  
  ainv=input.read<double>("aInv");
  
  //////////////////////////////////////////////////
  
  //set spatial sizes
  V=L[0]*pow(Ls,NDIM-1);
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  return *this;
}

perens_t& perens_t::allocate()
{
  for(auto &task : get_Zq_tasks(*this))
    task.out->resize(im_r_ilinmom_ind.max());
  
  for(auto &task : concat(get_Zbil_tasks(*this),get_pr_bil_tasks(*this)))
    task.out->resize(im_r_im_r_iZbil_ibilmom_ind.max());
  
  for(auto &task : concat(get_Zmeslep_tasks(*this),get_pr_meslep_tasks(*this)))
    task.out->resize(im_r_im_r_iop_iproj_imeslepmom_ind.max());

  return *this;
}

perens_t& perens_t::set_indices()
{
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  r_r_iZbil_ind.set_ranges({{"r",nr},{"r",nr},{"Zbil",nZbil}});
  im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  im_r_ijackp1_ind=im_r_ind*index_t({{"ijack",njacks+1}});
  
  im_r_im_r_iZbil_ind=im_r_ind*im_r_ind*index_t({{"iZbil",nZbil}});;
  r_ilinmom_ind.set_ranges({{"r",nr},{"linmom",linmoms.size()}});
  im_r_ilinmom_ind.set_ranges({{"m",nm},{"r",nr},{"linmom",linmoms.size()}});
  iZbil_ibilmom_ind.set_ranges({{"iZbil",nZbil},{"bilmom",bilmoms.size()}});
  im_r_im_r_iZbil_ibilmom_ind=im_r_im_r_iZbil_ind*index_t({{"bilmom",bilmoms.size()}});
  
  im_r_im_r_ilistGl_ipGl_ind=im_r_ind*im_r_ind*index_t({{"listGl",meslep::listGl.size()},{"ipGl",nGamma}});
  im_r_im_r_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"iop",nZbil},{"iproj",nZbil}});
  im_r_im_r_iop_iproj_imeslepmom_ind=im_r_im_r_iop_iproj_ind*index_t({{"imeslepmom",meslepmoms().size()}});
  im_r_im_r_iop_ilistpGl_ind.set_ranges({{"m_fw",nm},{"r_fw",nr},{"m_bw",nm},{"r_bw",nr},{"iop",meslep::nZop},{"listpGl",meslep::listpGl.size()}});
  iGl_ipGl_iclust_ind.set_ranges({{"iGamma",nGamma},{"ipGl",nGamma},{"iclust",njacks}});
  iop_ipGl_iclust_ind.set_ranges({{"iop",nZbil},{"ipGl",nGamma},{"iclust",njacks}});
  iop_iproj_iclust_ind.set_ranges({{"iop",nZbil},{"iproj",nZbil},{"iclust",njacks}});
  
  return *this;
}
