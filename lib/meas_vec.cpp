#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif
#include <meas_vec.hpp>

djvec_t read_conf_set_t(const string &template_path,vector<size_t> &id_list,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity)
{
  check_njacks_init();
  
  //get the list of existing files
  size_t clust_size=trim_to_njacks_multiple(id_list);
  cout<<"Clust size: "<<clust_size<<endl;
  
  //open all files
  vector<obs_file_t> files;
  for(auto &id : id_list) files.push_back(obs_file_t(combine(template_path.c_str(),id),ntot_col,cols));
  
  //measure the length of the first file
  size_t length;
  if(files.size()==0) length=0;
  else length=files[0].length(nlines);
  
  if(verbosity) cout<<"Total length (in multiple of nlines="<<nlines<<"): "<<length<<endl;
  
  //check that the length is a multiple of ncols*nlines
  size_t block_nentr=nlines*cols.size();
  size_t nblocks=length/block_nentr;
  if(length!=nblocks*block_nentr) CRASH("Total length is not a multiple of ncols*nlines=%zu",block_nentr);
  
  //! output
  cout<<"Allocating data"<<endl;
  djvec_t data(length);
  
  cout<<"Starting to read"<<endl;
#pragma omp parallel for
  for(size_t ijack=0;ijack<njacks;ijack++)
    {
      const size_t beg_file=ijack*clust_size;
      const size_t end_file=(ijack+1)*clust_size;
      
      if(verbosity)
	printf("Block of ijack: %zu, reading from file %zu to %zu\n",ijack,beg_file,end_file);
      
      for(size_t ifile=beg_file;ifile<end_file;ifile++)
	{
#ifdef USE_OMP
	  printf("Thread %d/%d reading file %zu/%zu\n",omp_get_thread_num(),omp_get_num_threads(),ifile,files.size());
#else
	  printf("Reading file %zu/%zu\n",ifile,files.size());
#endif
	  
	  //read all blocks
	  for(size_t iblock=0;iblock<nblocks;iblock++)
	    {
	      vector<double> temp=files[ifile].read(nlines);
	      if(temp.size()!=block_nentr) CRASH("Error reading file %zu, iblock %zu",ifile,iblock);
	      
	      //copy
	      for(size_t ientr=0;ientr<block_nentr;ientr++) data[ientr+block_nentr*iblock][ijack]+=temp[ientr];
	    }
	}
    }
  if(verbosity) cout<<"Finished reading"<<endl;
  
  //clusterize each entry
  {
    auto start=take_time();
    for(auto &d : data) d.clusterize(clust_size);
    if(verbosity) cout<<elapsed_time(start)<<" to clusterize"<<endl;
  }
  
  return data;
}

djvec_t read_conf_set_t(const string &template_path,const range_t &range,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity)
{
  //get the list of existing files
  vector<size_t> id_list=get_existing_paths_in_range(template_path,range,verbosity);
  
  return read_conf_set_t(template_path,id_list,ntot_col,cols,nlines,verbosity);
}
