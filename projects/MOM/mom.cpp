#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <corrections.hpp>
#include <geometry.hpp>
#include <types.hpp>

#include <Zbil.hpp>
#include <Zq.hpp>
#include <Zq_sig1.hpp>

//! compute the inverse of the passed propagator
vjprop_t get_prop_inv(const vjprop_t &jprop)
{
  vjprop_t jprop_inv(jprop.size());
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      put_into_jackknife(jprop_inv[imom],get_from_jackknife(jprop[imom],ijack).inverse(),ijack);
  
  return jprop_inv;
}

//! read a propagator from a given file
vprop_t read_prop(const string &path)
{
  vprop_t prop(imoms.size());
  
  //! source file
  raw_file_t file(path,"r");
  
  for(size_t is_so=0;is_so<NSPIN;is_so++)
    for(size_t ic_so=0;ic_so<NCOL;ic_so++)
      for(size_t imom=0;imom<imoms.size();imom++)
	for(size_t is_si=0;is_si<NSPIN;is_si++)
	  for(size_t ic_si=0;ic_si<NCOL;ic_si++)
	    file.bin_read(prop[imom](isc(is_si,ic_si),isc(is_so,ic_so)));
  
  return prop;
}

//! add the prop of a given conf on the jackknife
void build_jackknifed_prop(vjprop_t &jprop,const vprop_t &prop,size_t ijack)
{
  for(size_t imom=0;imom<imoms.size();imom++)
    add_to_cluster(jprop[imom],prop[imom],ijack);
}

//! compute all 16 vertices, with a given pair of propagators
void build_jackknifed_verts(vector<jverts_t> &jverts,const vprop_t &prop1,const vprop_t &prop2,size_t ijack)
{
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t iG=0;iG<nGamma;iG++)
      add_to_cluster(jverts[imom][iG],prop1[imom]*Gamma[iG]*Gamma[5]*prop2[imom].adjoint()*Gamma[5],ijack);
}

//! clusterize the propagator
void clusterize_prop(vjprop_t &jprop,size_t clust_size=1)
{
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    clusterize(jprop[imom],clust_size);
}

//! clusterize the vertexes
void clusterize_verts(vector<jverts_t> &jverts,size_t clust_size=1)
{
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t iG=0;iG<nGamma;iG++)
      clusterize(jverts[imom][iG],clust_size);
}

//! write a given Z
void write_Z(const string &name,const djvec_t &Z,const vector<double> &pt2)
{
  grace_file_t outf("plots/"+name+".xmg");
  outf<<fixed;
  outf.precision(8);
  outf.write_vec_ave_err(pt2,Z.ave_err());
}

int main(int narg,char **arg)
{
  //read input file
  string input_path="input.txt";
  if(narg>=2) input_path=arg[1];
  
  //open input
  raw_file_t input(input_path,"r");
  
  size_t Ls=input.read<size_t>("L"); //!< lattice spatial size
  L[0]=input.read<size_t>("T");
  
  const string act_str=input.read<string>("Action"); //!< action name
  auto act_key=gaz_decr.find(act_str); //!< key in the map of act
  if(act_key==gaz_decr.end()) CRASH("Unable to decript %s",act_str.c_str());
  gaz_t act=act_key->second; //!< action
  
  const double beta=input.read<double>("Beta"); //!< beta
  const double plaq=input.read<double>("Plaq"); //!< plaquette
  
  const string mom_list_path=input.read<string>("MomList"); //!< list of momenta
  const size_t ext_njacks=input.read<size_t>("NJacks"); //!< number of jacks
  
  //! conf range
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  //! template path
  string template_path=input.read<string>("TemplatePath");
  
  //////////////////////////////////////////////////
  
  //set the number of jackknives
  set_njacks(ext_njacks);
  
  double g2=6.0/beta; //!< coupling
  double g2tilde=g2/plaq; //!< boosted coupling
  
  //set the coefficients
  set_pr_bil_a2(act);
  
  //set spatial sizes
  V=L[0]*pow(Ls,NDIM-1);
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  //initialize momenta
  set_list_of_moms(mom_list_path);
  set_class_of_equiv_moms();
  set_filtered_moms();
  //list_all_smom_pairs();
  
  vector<size_t> conf_list=get_existing_paths_in_range(template_path,conf_range); //!< list of existing confs
  size_t clust_size=trim_to_njacks_multiple(conf_list,true); //!< cluster size
  vjprop_t jprop(imoms.size()); //!< jackkniffed propagator
  vector<jverts_t> jverts(imoms.size()); //!< jackkniffed vertex
  
#pragma omp parallel for
  for(size_t ijack=0;ijack<njacks;ijack++)
    for(size_t iconf=ijack*clust_size;iconf<(ijack+1)*clust_size;iconf++)
      {
	string path=combine(template_path.c_str(),conf_list[iconf]);
#ifdef USE_OMP
	printf("Thread %d/%d reading file %zu/%zu\n",omp_get_thread_num()+1,omp_get_num_threads(),iconf+1,conf_list.size());
#else
	printf("Reading file %zu/%zu\n",iconf+1,conf_list.size());
#endif
	
	vprop_t prop=read_prop(path); //!< propagator on a given conf
	build_jackknifed_prop(jprop,prop,ijack);
	build_jackknifed_verts(jverts,prop,prop,ijack);
      }
  
  //clusterize
  clusterize_prop(jprop,clust_size);
  clusterize_verts(jverts,clust_size);
  
  vjprop_t jprop_inv=get_prop_inv(jprop); //!< inverse prop
  
  //compute Zq, Zq_sig1 and Zbil for all moms
  djvec_t Zq_allmoms=compute_Zq(jprop_inv);
  djvec_t Zq_sig1_allmoms=compute_Zq_sig1(jprop_inv);
  vector<djvec_t> pr_bil_allmoms=compute_proj_bil(jprop_inv,jverts,jprop_inv);
  
  //correct Z
  djvec_t Zq_allmoms_sub=Zq_allmoms;
  djvec_t Zq_sig1_allmoms_sub=Zq_sig1_allmoms;
  vector<djvec_t> pr_bil_allmoms_sub=pr_bil_allmoms;
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      imom_t mom=imoms[imom];
      Zq_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
      Zq_sig1_allmoms_sub[imom]-=g2tilde*sig1_a2(act,mom,L);
      for(size_t iZbil=0;iZbil<nZbil;iZbil++)
	pr_bil_allmoms_sub[iZbil][imom]-=g2tilde*pr_bil_a2(act,mom,L,iZbil);
    }
  
  //average equiv moms
  djvec_t Zq=average_equiv_moms(Zq_allmoms);
  djvec_t Zq_sub=average_equiv_moms(Zq_allmoms_sub);
  djvec_t Zq_sig1=average_equiv_moms(Zq_sig1_allmoms);
  djvec_t Zq_sig1_sub=average_equiv_moms(Zq_sig1_allmoms_sub);
  vector<djvec_t> Zbil(nZbil);
  vector<djvec_t> Zbil_sub(nZbil);
  for(size_t iZbil=0;iZbil<nZbil;iZbil++)
    {
      Zbil[iZbil]=average_equiv_moms(Zq_allmoms/pr_bil_allmoms[iZbil]);
      Zbil_sub[iZbil]=average_equiv_moms(Zq_allmoms_sub/pr_bil_allmoms_sub[iZbil]);
    }
  
  //filter moms
  djvec_t Zq_sub_filt=get_filtered_moms(Zq_sub);
  djvec_t Zq_sig1_sub_filt=get_filtered_moms(Zq_sig1_sub);
  vector<djvec_t> Zbil_sub_filt(nZbil);
  for(size_t iZbil=0;iZbil<nZbil;iZbil++)
    Zbil_sub_filt[iZbil]=get_filtered_moms(Zbil_sub[iZbil]);
  
  write_Z("Zq",Zq,get_indep_pt2());
  write_Z("Zq_sub",Zq_sub,get_indep_pt2());
  write_Z("Zq_sig1",Zq_sig1,get_indep_pt2());
  write_Z("Zq_sig1_sub",Zq_sig1_sub,get_indep_pt2());
  write_Z("ZS",Zbil[iZS],get_indep_pt2());
  write_Z("ZA",Zbil[iZA],get_indep_pt2());
  write_Z("ZP",Zbil[iZP],get_indep_pt2());
  write_Z("ZV",Zbil[iZV],get_indep_pt2());
  write_Z("ZT",Zbil[iZT],get_indep_pt2());
  
  write_Z("Zq_sub_filt",Zq_sub_filt,get_filtered_pt2());
  write_Z("Zq_sig1_sub_filt",Zq_sig1_sub_filt,get_filtered_pt2());
  write_Z("ZS_sub_filt",Zbil_sub_filt[iZS],get_filtered_pt2());
  write_Z("ZA_sub_filt",Zbil_sub_filt[iZA],get_filtered_pt2());
  write_Z("ZP_sub_filt",Zbil_sub_filt[iZP],get_filtered_pt2());
  write_Z("ZV_sub_filt",Zbil_sub_filt[iZV],get_filtered_pt2());
  write_Z("ZT_sub_filt",Zbil_sub_filt[iZT],get_filtered_pt2());
  
  return 0;
}
