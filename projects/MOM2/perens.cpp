#define EXTERN_PERENS
 #include <MOM2/perens.hpp>

void perens_t::read_pars(const string &name)
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
}
