#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <geometry.hpp>

using vprop_t=vector<prop_t>;
using jverts_t=array<jprop_t,nGamma>;
using vjprop_t=vector<jprop_t>;

int main(int narg,char **arg)
{
  //read input file
  string input_path="input.txt";
  if(narg>=2) input_path=arg[1];
  
  //open input
  raw_file_t input(input_path,"r");
  
  //! lattice size
  size_t Ls=input.read<size_t>("L");
  L[0]=input.read<size_t>("T");
  
  //! list of momenta
  string mom_list_path=input.read<string>("MomList");
  
  //! number of jacks
  const size_t ext_njacks=input.read<size_t>("NJacks");
  
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
  
  //set spatial sizes
  V=L[0]*pow(Ls,NDIM-1);
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  //initialize momenta
  get_list_of_moms(mom_list_path);
  get_class_of_equiv_moms();
  //list_all_smom_pairs();
  
  //! list of existing confs
  vector<size_t> conf_list=get_existing_paths_in_range(template_path,conf_range);
  size_t clust_size=trim_to_njacks_multiple(conf_list,true);
  
  //! jackkniffed propagator
  vjprop_t jprop(imoms.size());
  
  //! jackkniffed vertex
  vector<jverts_t> jverts(imoms.size());

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
	
	//! source file
	raw_file_t file(path,"r");
	
	//! propagator on a given conf
	vprop_t prop(imoms.size());
	for(size_t is_so=0;is_so<NSPIN;is_so++)
	  for(size_t ic_so=0;ic_so<NCOL;ic_so++)
	    for(size_t imom=0;imom<imoms.size();imom++)
	      for(size_t is_si=0;is_si<NSPIN;is_si++)
		for(size_t ic_si=0;ic_si<NCOL;ic_si++)
		  file.bin_read(prop[imom](isc(is_si,ic_si),isc(is_so,ic_so)));
	
	for(size_t imom=0;imom<imoms.size();imom++)
	  {
	    // build the jackkniffed propagator
	    put_into_cluster(jprop[imom],prop[imom],ijack);
	    
	    //! compute all 16 vertices
	    for(size_t iG=0;iG<nGamma;iG++)
	      put_into_cluster(jverts[imom][iG],prop[imom]*Gamma[iG]*Gamma[5]*prop[imom].adjoint()*Gamma[5],ijack);
	  }
      }
  
  //////////////////////////////////////////// clusterize .///////////////////////////////////////
  
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    {
      clusterize(jprop[imom]);
      for(size_t iG=0;iG<nGamma;iG++) clusterize(jverts[imom][iG],clust_size);
    }
  
  //////////////////////////////// compute the inverse propagator ////////////////////////////////
  
  //! inverse prop
  vjprop_t jprop_inv(jprop.size());
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++)
    for(size_t ijack=0;ijack<=njacks;ijack++)
      put_into_jackknife(jprop_inv[imom],get_from_jackknife(jprop[imom],ijack).inverse(),ijack);
  
  //////////////////////////////////// compute Z by amputating //////////////////////////////////
  
  //! Z of quark field
  djvec_t Zq(equiv_imoms.size());
  
  //! sixteen gamma Z
  vector<djvec_t> Z(nGamma,djvec_t(equiv_imoms.size()));
  
#pragma omp parallel for
  for(size_t ind_mom=0;ind_mom<equiv_imoms.size();ind_mom++)
    {
#ifdef USE_OMP
      printf("Thread %d/%d analyzing mom %zu/%zu\n",omp_get_thread_num()+1,omp_get_num_threads(),ind_mom+1,equiv_imoms.size());
#else
      printf("Analyzing mom %zu/%zu\n",ind_mom+1,equiv_imoms.size());
#endif

      //reset Z
      Zq[ind_mom]=0;
      for(size_t iG=0;iG<nGamma;iG++) Z[iG][ind_mom]=0.0;
      
      //loop on equivalent moms
      auto &imom_class=equiv_imoms[ind_mom];
      double pt2=imoms[imom_class.first].p(L).tilde().norm2();
      for(size_t imom : imom_class.second)
	{
	  p_t ptilde=imoms[imom].p(L).tilde();
	  prop_t pslash=slash(ptilde);
	  
	  for(size_t ijack=0;ijack<=njacks;ijack++)
	    {
	      prop_t prop_inv=get_from_jackknife(jprop_inv[imom],ijack);
	      Zq[ind_mom][ijack]+=(prop_inv*pslash).trace().imag();
	      
	      for(size_t iG=0;iG<nGamma;iG++)
		{
		  prop_t vert=get_from_jackknife(jverts[imom][iG],ijack) ;
		  prop_t amp_vert=prop_inv*vert*Gamma[5]*prop_inv.adjoint()*Gamma[5];
		  Z[iG][ind_mom][ijack]+=(amp_vert*Gamma[iG]).trace().real();
		}
	    }
	}
      
      //normalize Z
      Zq[ind_mom]/=12*pt2*V*imom_class.second.size();
      for(size_t iG=0;iG<nGamma;iG++)
	Z[iG][ind_mom]/=12*imom_class.second.size();
    }
  
  //! normalize
  djvec_t ZS=Zq/Z[0];
  djvec_t ZA=4.0*Zq/(Z[1]+Z[2]+Z[3]+Z[4]);
  djvec_t ZP=Zq/Z[5];
  djvec_t ZV=4.0*Zq/(Z[6]+Z[7]+Z[8]+Z[9]);
  djvec_t ZT=6.0*Zq/(Z[10]+Z[11]+Z[12]+Z[13]+Z[14]+Z[15]);
  for(auto &p : vector<pair<djvec_t,string>>{{Zq,"Zq"},{ZS,"ZS"},{ZA,"ZA"},{ZP,"ZP"},{ZV,"ZV"},{ZT,"ZT"}})
    grace_file_t("plots/"+p.second+".xmg").write_vec_ave_err(get_indep_pt2(),p.first.ave_err());
  
  return 0;
}
