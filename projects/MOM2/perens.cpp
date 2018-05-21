#define EXTERN_PERENS
 #include <MOM2/perens.hpp>

#include <MOM2/sigma.hpp>
#include <MOM2/pr_meslep.hpp>

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
  
  im_sea=input.read<int>("ImSea");
  if(pars::chir_extr_method==chir_extr::MQUARK and (im_sea<0 or im_sea>=(int)nm))
    CRASH("When chiral extrapolation is done in terms of quark mass, im_sea must be in the range  [0,%zu), value %d is invalid",nm,im_sea);
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  
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
  
  set_comp_list_of_moms(dir_path+"/"+mom_list_path,pars::filter_thresh);
  
  return *this;
}

void perens_t::write_pars(const string &name) const
{
  raw_file_t output(name+"/pars.txt","w");
  
  output.write(L[1],"L");
  output.write(L[0],"T");
  
  output.write(beta,"Beta");
  output.write(plaq,"Plaq");
  output.write(nm,"Nm");
  for(size_t im=0;im<nm;im++) output.write(am[im]);
  output.write(nr,"Nr");
  output.write(im_sea,"ImSea");
  
  output.write("mom_list.txt","MomList");
  
  output.write("ConfRange");
  output.write(conf_range.start);
  output.write(conf_range.each);
  output.write(conf_range.end);
  
  output.write(prop_hadr_path,"PropHadrPath");
  output.write(prop_lep_path,"PropLepPath");
  
  output.write(nhits,"NHits");
  output.write(nhits_to_use,"NHitsToUse");
  
  output.write(tmin,"TMin");
  output.write(tmax,"TMax");
  
  output.write(tmax,"TMax");
  output.write(ainv,"ainv");
}

perens_t& perens_t::allocate()
{
  for(auto &Ztask : get_all_Ztasks())
    Ztask.out->resize(Ztask.ind.max());
  
  sigma.resize(im_r_ilinmom_isigmaproj_isigmains_ind.max());
  pr_bil.resize(im_r_im_r_bilins_ibil_ibilmom_ind.max());
  pr_meslep.resize(im_r_im_r_meslepins_iop_iproj_imeslepmom_ind.max());
  
  for(auto &t : {&deltam_cr,&deltam_tm})
    t->resize(im_r_ind.max());
  
  meson_mass.resize(im_im_ind.max());
  
  return *this;
}

perens_t& perens_t::set_indices()
{
  im_im_ind.set_ranges({{"m",nm},{"m",nm}});
  
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  im_r_im_r_bilins_igam_ind=im_r_ind*im_r_ind*index_t({{"bilins",pr_bil::nins},{"igamma",nGamma}});
  im_r_im_r_bilins_ibil_ibilmom_ind=im_r_ind*im_r_ind*index_t({{"bilins",pr_bil::nins},{"bil",nbil},{"bilmom",bilmoms.size()}});
  
  im_r_im_r_meslepins_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"meslepins",pr_meslep::nins},{"iop",nbil},{"iproj",nbil}});;
  im_r_im_r_meslepins_iop_iproj_imeslepmom_ind=im_r_im_r_meslepins_iop_iproj_ind*index_t({{"imeslepmom",meslepmoms().size()}});
  
  r_r_ibil_ind.set_ranges({{"r",nr},{"r",nr},{"bil",nbil}});
  im_r_ijack_ind=im_r_ind*index_t({{"ijack",njacks}});
  im_r_ijackp1_ind=im_r_ind*index_t({{"ijack",njacks+1}});
  
  im_r_ilinmom_isigmaproj_isigmains_ind.set_ranges({{"m",nm},{"r",nr},{"linmom",linmoms.size()},{"proj",sigma::nproj},{"ins",sigma::nins}});
  
  im_r_im_r_ibil_ind=im_r_ind*im_r_ind*index_t({{"ibil",nbil}});;
  r_ilinmom_ind.set_ranges({{"r",nr},{"linmom",linmoms.size()}});
  im_r_ilinmom_ind.set_ranges({{"m",nm},{"r",nr},{"linmom",linmoms.size()}});
  ibil_ibilmom_ind.set_ranges({{"ibil",nbil},{"bilmom",bilmoms.size()}});
  im_r_im_r_ibil_ibilmom_ind=im_r_im_r_ibil_ind*index_t({{"bilmom",bilmoms.size()}});
  
  im_r_im_r_ilistGl_ipGl_ind=im_r_ind*im_r_ind*index_t({{"listGl",meslep::listGl.size()},{"ipGl",nGamma}});
  im_r_im_r_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"iop",nbil},{"iproj",nbil}});
  im_r_im_r_iop_iproj_imeslepmom_ind=im_r_im_r_iop_iproj_ind*index_t({{"imeslepmom",meslepmoms().size()}});
  im_r_im_r_iop_ilistpGl_ind.set_ranges({{"m_fw",nm},{"r_fw",nr},{"m_bw",nm},{"r_bw",nr},{"iop",meslep::nZop},{"listpGl",meslep::listpGl.size()}});
  iGl_ipGl_iclust_ind.set_ranges({{"iGamma",nGamma},{"ipGl",nGamma},{"iclust",njacks}});
  iop_ipGl_iclust_ind.set_ranges({{"iop",nbil},{"ipGl",nGamma},{"iclust",njacks}});
  iop_iproj_iclust_ind.set_ranges({{"iop",nbil},{"iproj",nbil},{"iclust",njacks}});
  
  ilistGl_ilistpGl_iclust_ind=index_t({{"listGl",meslep::listGl.size()},{"pGl",meslep::listpGl.size()},{"iclust",njacks}});
  ilistGl_ilistpGl_lins_iclust_ind=index_t({{"listGl",meslep::listGl.size()},{"pGl",meslep::listpGl.size()},{"lins",lprop::nins},{"iclust",njacks}});
  im_r_im_r_iop_ilistpGl_meslepins_ind.set_ranges({{"m_fw",nm},{"r_fw",nr},{"m_bw",nm},{"r_bw",nr},{"iop",meslep::nZop},{"listpGl",meslep::listpGl.size()},{"meslepins",pr_meslep::nins}});
  
  return *this;
}

perens_t& perens_t::set_pars_for_scratch()
{
  using namespace reno_scheme;
  
  switch(pars::scheme)
    {
    case RI_MOM:
      set_ri_mom_moms();
      break;
    case SMOM:
      set_smom_moms();
      break;
    }
  
  //write all linear momenta
  const string linpath=dir_path+"/linmoms.txt";
  ofstream linmom_file(linpath);
  if(not linmom_file.good())
    CRASH("Failed to open %s",linpath.c_str());
  for(const auto &m : linmoms)
    linmom_file<<m[0]<<"\t=\t"<<all_moms[m[0]]<<endl;
  
  //write all bilinear combo
  const string bilpath=dir_path+"/bilmoms.txt";
  ofstream bilmom_file(bilpath);
  if(not bilmom_file.good())
    CRASH("Failed to open %s",bilpath.c_str());
  for(const auto &m : bilmoms)
    {
      const size_t mom=m[0];
      const size_t ilinmom1=m[1];
      const size_t ilinmom2=m[2];
      const size_t mom1=linmoms[ilinmom1][0];
      const size_t mom2=linmoms[ilinmom2][0];
      bilmom_file<<mom<<all_moms[mom]<<"\t"
		 <<mom1<<all_moms[mom1]<<"\t"
		 <<mom2<<all_moms[mom2]<<endl;
    }
  
  return *this;
}

void perens_t::prepare_list_of_confs()
{
  //create the list of all confs available
  string test_path=dir_path+"/"+prop_hadr_path+"/%04zu/fft_S_M0_R0_0";
  if(nhits>1) test_path+="_hit_0";
  
  conf_list=get_existing_paths_in_range(test_path,conf_range);
  if(conf_list.size()==0) CRASH("list of configurations is empty! check %s ",test_path.c_str());
  clust_size=trim_to_njacks_multiple(conf_list,true);
  
  conf_ind.set_ranges({{"ijack",njacks},{"i_in_clust",clust_size}});
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  i_in_clust_ihit_ind.set_ranges({{"i_in_clust",clust_size},{"ihit",nhits_to_use}});
  im_r_iconf_ihit_iqins_ind=im_r_ind*index_t({{"conf",conf_list.size()},{"hit",nhits},{"qins",qprop::nins}});
  im_r_iqins_ijack_ind=im_r_ind*index_t({{"qins",qprop::nins},{"ijack",njacks}});
  im_r_ijqins_ijack_ind=im_r_ind*index_t({{"jqins",jqprop::nins},{"ijack",njacks}});
  im_r_ijqins_ind=im_r_ind*index_t({{"jqins",jqprop::nins}});
  ilins_ijack_ind=index_t({{"lins",lprop::nins},{"ijack",njacks}});
  iconf_ihit_ilins_ind=index_t({{"conf",conf_list.size()},{"hit",nhits},{"lins",lprop::nins}});
}

perens_t& perens_t::compute_ingredients(const string& ingredients_path)
{
  prepare_list_of_confs();
  
  using namespace reno_scheme;
  
  switch(pars::scheme)
    {
    case RI_MOM:
      mom();
      break;
    case SMOM:
      CRASH("Here or deeper branching?");
      //smom();
      break;
    }
  bin_write_ingredients(ingredients_path);
  
  return *this;
}

perens_t& perens_t::read_or_compute_ingredients()
{
  const string ingredients_path=dir_path+"/ingredients.dat";
  
  //if ingredients exists read it, otherwise compute it
  if(file_exists(ingredients_path))
    {
      cout<<"File "<<ingredients_path<<" found, opening"<<endl;
      bin_read_ingredients(ingredients_path);
    }
  else
    {
      cout<<"File "<<ingredients_path<<" not found, computing"<<endl;
      compute_ingredients(ingredients_path);
    }
  
  return *this;
}

void perens_t::bin_read_ingredients(raw_file_t &file)
{
  for(auto t : {&sigma,&pr_bil,&pr_meslep})
    t->bin_read(file);
}

void perens_t::bin_write_ingredients(raw_file_t &file)
{
  for(auto t : {&sigma,&pr_bil,&pr_meslep})
    t->bin_write(file);
}

perens_t perens_t::average_r() const
{
  perens_t out=*this;
  
  if(nr==1)
      cout<<"Skipping r average"<<endl;
  else
    {
      out.nr=1;
      out.linmoms=linmoms;
      out.bilmoms=bilmoms;
      
      out.set_indices();
      out.allocate();
      
      average_r_sigma(out);
      average_r_pr_bil(out);
      average_r_pr_meslep(out);
      average_r_deltam(out);
      
      flush_unused_memory();
  }
  
  return out;
}

perens_t perens_t::val_chir_extrap() const
{
  perens_t out=*this;
  
  if(nm>1)
    {
      out.nm=1;
      out.am={0.0};
      out.meson_mass_sea=meson_mass_sea;
      out.meson_mass=0.0;
      
      out.set_indices();
      out.allocate();
      
      val_chir_extrap_sigma(out);
      val_chir_extrap_pr_bil(out);
      val_chir_extrap_pr_meslep(out);
      
      val_chir_extrap_deltam(out);
    }
  else
    cout<<"Skipping Valence chiral extrapolation "<<dir_path<<endl;
  
  return out;
}

perens_t perens_t::average_equiv_momenta() const
{
  perens_t out=*this;
  
  //build out lin list
  const vector<vector<size_t>> equiv_linmom_combos=get_equiv_list(linmoms,dir_path+"/equiv_linmoms.txt");
  fill_output_equivalent_momenta(out.linmoms,equiv_linmom_combos,equiv_linmom_combos,linmoms);
  
  //build out bil combo
  const vector<vector<size_t>> equiv_bilmom_combos=get_equiv_list(bilmoms,dir_path+"/equiv_bilmoms.txt");
  fill_output_equivalent_momenta(out.bilmoms,equiv_linmom_combos,equiv_bilmom_combos,bilmoms);
  
  //build out bil combo
  const vector<vector<size_t>> equiv_meslepmom_combos=get_equiv_list(meslepmoms(),dir_path+"/equiv_meslepmoms.txt");
  //fill_output_equivalent_momenta(out.meslepmoms(),equiv_linmom_combos,equiv_meslepmom_combos,meslepmoms());
  
  out.set_indices();
  out.allocate();
  
  average_equiv_momenta_sigma(out,equiv_linmom_combos);
  if(pars::compute_bilinears) average_equiv_momenta_pr_bil(out,equiv_bilmom_combos);
  if(pars::compute_meslep)    average_equiv_momenta_pr_meslep(out,equiv_meslepmom_combos);
  
  flush_unused_memory();
  
  return out;
}

void perens_t::print_discr()
{
  ofstream out(dir_path+"/p2_p4.txt");
  for(size_t imom=0;imom<linmoms.size();imom++)
    {
      const p_t p=all_moms[linmoms[imom][0]].p(L);
      const double p2tilde=p.tilde().norm2();
      out<<p2tilde<<" "<<p.norm2()<<" "<<p.norm4()<<endl;
    }
}

perens_t perens_t::interpolate_to_p2ref()
{
  perens_t out=*this;
  
  //reset momenta ranges
  out.all_moms=vector<imom_t>{1};
  out.linmoms=vector<array<size_t,1>>{{0}};
  out.bilmoms=vector<array<size_t,3>>{{0,0,0}};
  
  //set indices
  out.set_indices();
  out.allocate();
  
  interpolate_Zq_to_p2ref(out);
  if(pars::compute_bilinears) interpolate_Zbil_to_p2ref(out);
  if(pars::compute_meslep) interpolate_Zmeslep_to_p2ref(out);
  
  return out;
}


perens_t perens_t::write_checkpoint() const
{
  write_pars(dir_path+"/check");
  write_comp_list_of_moms(dir_path+"/check/momlist.txt");
  
  return *this;
}
