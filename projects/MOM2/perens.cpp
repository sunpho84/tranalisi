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

void perens_t::write_pars(const string &path) const
{
  ofstream output(path);
  if(not output.good()) CRASH("Unable to open %s",path.c_str());
  
  output<<"L "<<L[1]<<endl;
  output<<"T "<<L[0]<<endl;
  
  output<<"Beta "<<beta<<endl;
  output<<"Plaq "<<plaq<<endl;
  output<<"Nm "<<nm<<"   ";
  for(size_t im=0;im<nm;im++) output<<am[im]<<" ";
  output<<endl;
  
  output<<"NR "<<nr<<endl;
  output<<"ImSea "<<im_sea<<endl;
  
  output<<"MomList "<<"mom_list.txt"<<endl;
  
  output<<"ConfRange "<<" "<<conf_range.start<<" "<<conf_range.each<<" "<<conf_range.end<<endl;
  
  output<<"PropHadrPath "<<prop_hadr_path<<endl;
  output<<"PropLepPath "<<prop_lep_path<<endl;
  
  output<<"NHits "<<nhits<<endl;
  output<<"NHitsToUse "<<nhits_to_use<<endl;
  
  output<<"TMin "<<tmin<<endl;
  output<<"TMax "<<tmax<<endl;
  
  output<<"ainv "<<ainv<<endl;
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
  if(pars::use_QED) meson_mass_QED.resize(im_im_ind.max());
  
  return *this;
}

perens_t& perens_t::set_indices()
{
  im_im_ind.set_ranges({{"m",nm},{"m",nm}});
  
  im_r_ind.set_ranges({{"m",nm},{"r",nr}});
  im_r_im_r_igam_ind=im_r_ind*im_r_ind*index_t({{"igamma",nGamma}});
  im_r_im_r_bilins_igam_ind=im_r_ind*im_r_ind*index_t({{"bilins",pr_bil::nins},{"igamma",nGamma}});
  im_r_im_r_bilins_ibil_ibilmom_ind=im_r_ind*im_r_ind*index_t({{"bilins",pr_bil::nins},{"bil",nbil},{"bilmom",bilmoms.size()}});
  
  im_r_im_r_meslepins_iop_iproj_ind=im_r_ind*im_r_ind*index_t({{"meslepins",pr_meslep::nins},{"iop",meslep::nZop},{"iproj",meslep::nZop}});;
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

perens_t& perens_t::compute_ingredients()
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
  bin_write_ingredients();
  
  return *this;
}

perens_t& perens_t::read_or_compute_ingredients()
{
  //if ingredients exists read it, otherwise compute it
  if(file_exists(ingredients_path()))
    {
      cout<<"File "<<ingredients_path()<<" found, opening"<<endl;
      bin_read_ingredients();
    }
  else
    {
      cout<<"File "<<ingredients_path()<<" not found, computing"<<endl;
      compute_ingredients();
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
      cout<<meson_mass<<endl;
      cout<<out.meson_mass<<endl;
      
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

perens_t perens_t::interpolate_or_extrapolate_to_single_p2_preamble() const
{
  
  perens_t out=*this;
  
  //reset momenta ranges
  out.all_moms=vector<imom_t>{1};
  out.linmoms=vector<array<size_t,1>>{{0}};
  out.bilmoms=vector<array<size_t,3>>{{0,0,0}};
  
  //set indices
  out.set_indices();
  out.allocate();
  
  return out;
}

perens_t perens_t::interpolate_to_p2ref() const
{
  perens_t out=interpolate_or_extrapolate_to_single_p2_preamble();
  
  interpolate_Zq_to_p2ref(out);
  if(pars::compute_bilinears) interpolate_Zbil_to_p2ref(out);
  if(pars::compute_meslep) interpolate_Zmeslep_to_p2ref(out);
  
  return out;
}


perens_t perens_t::extrapolate_to_0_p2() const
{
  perens_t out=interpolate_or_extrapolate_to_single_p2_preamble();
  
  for(auto t : out.get_Zq_tasks({this}))
    extrapolate_to_0_p2_internal(linmoms,t);
  
  for(auto t : out.get_Zbil_tasks({this}))
    extrapolate_to_0_p2_internal(bilmoms,t);
  
  for(auto t : out.get_Zmeslep_tasks({this}))
    extrapolate_to_0_p2_internal(meslepmoms(),t);
  
  return out;
}

perens_t perens_t::assemble_QED_greenfunctions() const
{
  perens_t out=*this;
  
  out.assemble_sigma_QED_greenfunctions();
  if(pars::compute_bilinears) out.assemble_pr_bil_QED_greenfunctions();
  if(pars::compute_meslep) out.assemble_pr_meslep_QED_greenfunctions();
  
  return out;
}

perens_t perens_t::evolve() const
{
  perens_t out=*this;
  
  evolve_sigma(out);
  evolve_pr_bil(out);
  evolve_pr_meslep(out);
  
  return out;
}

perens_t perens_t::subtract_Oa2() const
{
  perens_t out=*this;
  
  out.subtract_Oa2_sigma();
  if(pars::compute_bilinears) out.subtract_Oa2_pr_bil();
  //subtract_Oa2_pr_meslep();
  
  return out;
}

perens_t perens_t::match_to_W_reg() const
{
  perens_t out=*this;
  
  match_to_W_reg_Zmeslep(out);
  
  return out;
}

void perens_t::print_Z(ofstream& file)
{
  using namespace pars;
  
  if(nm!=1) CRASH("Can print only 1 m, %zu present",nm);
  if(nr!=1) CRASH("Can print only 1 r, %zu present",nr);
  
  /////////////////////////////////////////////////////////////////
  file<<"/////////////////////////////////////////////////////////////////"<<endl;
  if(linmoms.size()!=1) CRASH("Can print only 1 mom, %zu present",linmoms.size());
  
  for(auto t : get_Zq_tasks())
    {
      file<<" "<<t.tag<<endl;
      file<<smart_print((*t.out)[0])<<endl;
      file<<endl;
    }
  
  /////////////////////////////////////////////////////////////////
  file<<"/////////////////////////////////////////////////////////////////"<<endl;
  if(bilmoms.size()!=1) CRASH("Can print only 1 mom, %zu present",bilmoms.size());
  
  for(auto t : get_Zbil_tasks())
    {
      file<<" "<<t.tag<<endl;
      for(size_t ibil=0;ibil<nbil;ibil++)
	file<<" "<<bil_tag[ibil]<<": "<<smart_print((*t.out)[ibil])<<endl;
      file<<endl;
    }
  
  /////////////////////////////////////////////////////////////////
  file<<"/////////////////////////////////////////////////////////////////"<<endl;
  if(bilmoms.size()!=1) CRASH("Can print only 1 mom, %zu present",meslepmoms().size());
  
  for(auto t : get_Zmeslep_tasks())
    {
      file<<" "<<t.tag<<endl;
      for(size_t iop=0;iop<nbil;iop++)
	{
	  for(size_t iproj=0;iproj<nbil;iproj++)
	    {
	      ave_err_t ae=(*t.out)[im_r_im_r_iop_iproj_imeslepmom_ind({0,0,0,0,iop,iproj,0})].ave_err();
	      if(isnan(ae.ave()))
		file<<"\t";
	      else
		file<<" "<<smart_print(ae);
	    }
	file<<endl;
	}
      file<<endl;
    }
}

perens_t perens_t::write_checkpoint()
{
  cout<<"Writing checkpoint for ens "<<dir_path<<endl;
  
  //take note of the old dir_path
  const string old_path=dir_path;
  dir_path=dir_path+"/check";
  
  //write everything down
  write_pars(dir_path+"/pars.txt");
  write_comp_list_of_moms(dir_path+"/mom_list.txt");
  bin_write_ingredients();
  bin_write_deltam();
  bin_write_meson_mass();
  bin_write_meson_mass_QED();
  bin_write_meson_mass_sea();
  
  //put back the old path
  dir_path=old_path;
  
  return *this;
}

void perens_t::make_Z_QED_absolute()
{
  if(not pars::use_QED)
    CRASH("Need to have been instructed to compute QED");
  
  Zq_QED_rel*=Zq;
  if(pars::compute_bilinears) Zbil_QED_rel*=Zbil;
  if(pars::compute_meslep) make_Zmeslep_QED_absolute();
}

void perens_t::plot_Z(const string &suffix)
{
  plot_sigma(suffix);
  plot_Zq(suffix);
  if(pars::compute_bilinears) plot_Zbil(suffix);
  if(pars::compute_meslep) plot_Zmeslep(suffix);
}
