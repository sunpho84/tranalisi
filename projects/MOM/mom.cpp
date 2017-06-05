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
  vprop_t prop(file_list.size(),vector<prop_t>(moms.size()));
 
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
	   for(size_t imom=0;imom<moms.size();imom++)
	     for(size_t is_si=0;is_si<NSPIN;is_si++)
	       for(size_t ic_si=0;ic_si<NCOL;ic_si++)
		 file.bin_read(prop[ifile][imom](isc(is_so,ic_so),isc(is_si,ic_si)));
   }
 
 return prop;
}

//! build the jackkniffed propagator
vjprop_t get_jprop(const vprop_t &prop,size_t clust_size)
{
  vjprop_t jprop(moms.size());
  for(size_t iconf=0;iconf<prop.size();iconf++)
#pragma omp parallel for
    for(size_t imom=0;imom<moms.size();imom++)
      put_into_cluster(jprop[imom],prop[iconf][imom],iconf/clust_size);
  for(auto &j : jprop) clusterize(j,clust_size);
  
  return jprop;
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
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  
  //initialize momenta
  get_list_of_moms(mom_list_path);
  get_class_of_equiv_moms();
  
  //! list of existing confs
  vector<size_t> conf_list=get_existing_paths_in_range(template_path,conf_range);
  size_t clust_size=trim_to_njacks_multiple(conf_list,true);
  
  //! one entry per file
  auto prop=read_prop(template_path,conf_list);
  cout<<"Finished reading"<<endl;
  
  //! jackkniffed propagator
  auto jprop=get_jprop(prop,clust_size);
  
  cout<<jprop[0]<<endl;
  
  // cout<<Gamma[5]<<endl;
  
  
  // //a+=b;
  
  // RowVectorXi adew;
  // //no no no fai jack l'indice più interno, punto
  // //cout<<a<<endl;
  // cout<<Gamma[4]<<endl;
  
  // //cout<<del<<endl;
  
  // //testing complex jack
  // jack_t<dcomplex> d,e;
  // d+=e;

  
  return 0;
}
