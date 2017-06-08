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
  
  //! Z
  enum{iZS,iZA,iZP,iZV,iZT};
  const size_t nZ=5,iZ_of_iG[nGamma]={0,1,1,1,1,2,3,3,3,3,4,4,4,4,4,4};
  const size_t Zdeg[nZ]={1,4,1,4,6};
  vector<djvec_t> Z(nZ,djvec_t(equiv_imoms.size()));
  
#pragma omp parallel for
  for(size_t ind_mom=0;ind_mom<equiv_imoms.size();ind_mom++)
    {
#ifdef USE_OMP
      printf("Thread %d/%d analyzing mom %zu/%zu\n",omp_get_thread_num()+1,omp_get_num_threads(),ind_mom+1,equiv_imoms.size());
#else
      printf("Analyzing mom %zu/%zu\n",ind_mom+1,equiv_imoms.size());
#endif

      //reset Z
      Zq[ind_mom]=0.0;
      for(auto &Zi : Z) Zi[ind_mom]=0.0;
      
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
	      double Zq_cl=(prop_inv*pslash).trace().imag()/(12.0*pt2*V);
	      Zq[ind_mom][ijack]+=Zq_cl;
	      
	      for(size_t iG=0;iG<nGamma;iG++)
		{
		  prop_t vert=get_from_jackknife(jverts[imom][iG],ijack) ;
		  prop_t amp_vert=prop_inv*vert*Gamma[5]*prop_inv.adjoint()*Gamma[5];
		  double pr=(amp_vert*Gamma[iG]).trace().real()/12.0;
		  double Z_cl=Zq_cl/pr;
		  
		  Z[iZ_of_iG[iG]][ind_mom][ijack]+=Z_cl;
		}
	    }
	}
      
      //normalize Z
      Zq[ind_mom]/=imom_class.second.size();
      for(size_t iZ=0;iZ<nZ;iZ++)
	Z[iZ][ind_mom]/=imom_class.second.size()*Zdeg[iZ];
    }
  
  //! normalize
  for(auto &p : vector<pair<djvec_t,string>>{{Zq,"Zq"},{Z[iZS],"ZS"},{Z[iZA],"ZA"},{Z[iZP],"ZP"},{Z[iZV],"ZV"},{Z[iZT],"ZT"}})
    {
      grace_file_t outf("plots/"+p.second+".xmg");
      outf<<fixed;
	    outf.precision(8);
	    outf.write_vec_ave_err(get_indep_pt2(),p.first.ave_err());
    }
	
	return 0;
}
