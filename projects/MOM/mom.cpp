#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <tranalisi.hpp>

#include <geometry.hpp>

using vprop_t=vector<vector<prop_t>>;
using vjprop_t=vector<jprop_t>;

//! read the propagator
vprop_t read_prop(const string &template_path,const vector<size_t> &file_list,bool verbosity=VERBOSE)
{
  auto start=take_time();
  
  vprop_t prop(file_list.size(),vector<prop_t>(imoms.size()));
  
#pragma omp parallel for
  for(size_t ifile=0;ifile<file_list.size();ifile++)
    {
      raw_file_t file(combine(template_path.c_str(),file_list[ifile]),"r");
      
      if(verbosity)
#ifdef USE_OMP
	printf("Thread %d/%d reading file %zu/%zu\n",omp_get_thread_num(),omp_get_num_threads(),ifile,file_list.size());
#else
        printf("Reading file %zu/%zu\n",ifile,files.size());
#endif
	
	for(size_t is_so=0;is_so<NSPIN;is_so++)
	  for(size_t ic_so=0;ic_so<NCOL;ic_so++)
	    for(size_t imom=0;imom<imoms.size();imom++)
	      for(size_t is_si=0;is_si<NSPIN;is_si++)
		for(size_t ic_si=0;ic_si<NCOL;ic_si++)
		  file.bin_read(prop[ifile][imom](isc(is_si,ic_si),isc(is_so,ic_so)));
    }
  
  cout<<"Finished reading in "<<elapsed_time(start)<<endl;
  
  return prop;
}

//! build the jackkniffed propagator
vjprop_t get_jprop(const vprop_t &prop,size_t clust_size)
{
  auto start=take_time();
  
  vjprop_t jprop(imoms.size());
  for(size_t iconf=0;iconf<prop.size();iconf++)
#pragma omp parallel for
    for(size_t imom=0;imom<imoms.size();imom++)
      put_into_cluster(jprop[imom],prop[iconf][imom],iconf/clust_size);
  for(auto &j : jprop) clusterize(j,clust_size);
  
  cout<<"Finished computing prop in "<<elapsed_time(start)<<endl;
  
  return jprop;
}

//! build the inverse
vjprop_t get_inv_jprop(const vjprop_t &jprop)
{
  auto start=take_time();
  
  vjprop_t jprop_inv(jprop.size());
#pragma omp parallel for
  for(size_t imom=0;imom<imoms.size();imom++) jprop_inv[imom]=jprop[imom].inverse();
  
  cout<<"Finished inverting the prop in "<<elapsed_time(start)<<endl;
  
  return jprop_inv;
}

//! compute Zq
djvec_t compute_Zq(const vjprop_t &jprop_inv)
{
  auto start=take_time();
  
  djvec_t out(equiv_imoms.size());
  
  //loop on equivalence class
#pragma omp parallel for
  for(size_t ind_mom=0;ind_mom<equiv_imoms.size();ind_mom++)
    {
      auto &imom_class=equiv_imoms[ind_mom];
      djack_t Zq=0;
      double pt2=imoms[imom_class.first].p(L).tilde().norm2();
      
      //loop on equivalent moms
      for(size_t imom : imom_class.second)
	{
	  p_t ptilde=imoms[imom].p(L).tilde();
	  prop_t pslash=slash(ptilde);
	  Zq+=(jprop_inv[imom]*pslash.cast<cdjack_t>()).trace().imag();
	}
      out[ind_mom]=Zq/(12*pt2*V*imom_class.second.size());
    }
  
  cout<<"Computed Zq in "<<elapsed_time(start)<<endl;
  
  return out;
}

//! compute projected amputated vertex on a given momentum
djvec_t compute_projected_amputated_vertex(const vprop_t &prop,const vjprop_t &prop_inv,const prop_t &G,const prop_t &Gproj,bool ri,size_t clust_size=1)
{
  auto start=take_time();
  
  djvec_t out(equiv_imoms.size());
  
  //loop on equivalence class
#pragma omp parallel for
  for(size_t ind_mom=0;ind_mom<equiv_imoms.size();ind_mom++)
    {
      auto &imom_class=equiv_imoms[ind_mom];
      out[ind_mom]=0;
      
      //loop on equivalent moms
      for(size_t imom : imom_class.second)
	{
	  jprop_t vert;
	  for(size_t iconf=0;iconf<prop.size();iconf++)
	    put_into_cluster(vert,prop[iconf][imom]*G*Gamma[5]*prop[iconf][imom].adjoint()*Gamma[5],iconf/clust_size);
	  clusterize(vert,clust_size);
	  
	  jprop_t amp_vert=prop_inv[imom]*vert*Gamma[5].cast<cdjack_t>()*prop_inv[imom].adjoint()*Gamma[5].cast<cdjack_t>();
	  cdjack_t pr=(amp_vert*Gproj.cast<cdjack_t>()).trace();
	  
	  out[ind_mom]+=get_re_or_im(pr,ri);
	}
      
      out[ind_mom]/=12*imom_class.second.size();
    }
  
  cout<<"Finished computing Z in "<<elapsed_time(start)<<endl;
  
  return out;
}

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
  
  //! one entry per file
  vprop_t prop=read_prop(template_path,conf_list);
  
  //! jackkniffed propagator
  vjprop_t jprop=get_jprop(prop,clust_size);
  
  //! inverse prop
  vjprop_t jprop_inv=get_inv_jprop(jprop);
  
  //! Zq
  djvec_t Zq=compute_Zq(jprop_inv);
  
  //! compute all 16 Z
  vector<djvec_t> Z(nGamma);
  for(size_t iG=0;iG<nGamma;iG++) Z[iG]=compute_projected_amputated_vertex(prop,jprop_inv,Gamma[iG],Gamma[iG],RE,clust_size);
  
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
