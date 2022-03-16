#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <index.hpp>
#include <meas_vec.hpp>

djvec_t read_conf_set_t(const vector<string> &file_paths,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity)
{
  //open all files
  vector<obs_file_t> files;
  for(auto &file_path : file_paths)
    files.push_back(obs_file_t(file_path,ntot_col,cols));
  
  const size_t nFiles=file_paths.size();
  
  //measure the length of the first file
  size_t length;
  if(nFiles==0) length=0;
  else length=files[0].length(nlines);
  
  if(verbosity==VERBOSE) cout<<"Total length (in multiple of nlines="<<nlines<<"): "<<length<<endl;
  
  const double clust_size=nFiles/(double)njacks;
  
  if(verbosity==VERBOSE) cout<<"Clust size: "<<clust_size<<endl;
  
  //check that the length is a multiple of ncols*nlines
  size_t block_nentr=nlines*cols.size();
  size_t nblocks=length/block_nentr;
  if(length!=nblocks*block_nentr) CRASH("Total length is not a multiple of ncols*nlines=%zu",block_nentr);
  
  //! output
  if(verbosity==VERBOSE) cout<<"Allocating data"<<endl;
  djvec_t data(length);
  
  if(verbosity==VERBOSE)cout<<"Starting to read"<<endl;
  
  index_t index({{"file",nFiles},{"block",nblocks},{"block_entr",block_nentr}});
  
  vector<double> rawData;
  
#pragma omp parallel for
  for(size_t ifile=0;ifile<nFiles;ifile++)
    {
      if(verbosity)
#ifdef USE_OMP
	printf("Thread %d/%d reading file %zu/%zu\n",omp_get_thread_num(),omp_get_num_threads(),ifile,files.size());
#else
      printf("Reading file %zu/%zu\n",ifile,files.size());
#endif
      
      //read all blocks
      for(size_t iblock=0;iblock<nblocks;iblock++)
	{
	  vector<double> temp=files[ifile].read(nlines,false);
	  if(temp.size()!=block_nentr) CRASH("Error reading file %zu, iblock %zu",ifile,iblock);
	  
	  //copy
	  for(size_t ientr=0;ientr<block_nentr;ientr++) rawData[index({ifile,iblock,ientr})]=temp[ientr];
	}
    }
    
  if(verbosity) cout<<"Finished reading"<<endl;
  
  jackknivesFill(nFiles,[&](const size_t& iConf,const size_t& iClust,const double& weight)
  {
    for(size_t iblock=0;iblock<nblocks;iblock++)
      for(size_t ientr=0;ientr<block_nentr;ientr++)
	data[ientr+block_nentr*iblock][iClust]+=rawData[index({iConf,iblock,ientr})];
  });
  
  //clusterize each entry
  {
    auto start=take_time();
    for(auto &d : data) d.clusterize(clust_size);
    if(verbosity) cout<<elapsed_time(start)<<" to clusterize"<<endl;
  }
  
  return data;
}

djvec_t read_conf_set_t(const string &template_path,vector<size_t> &id_list,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity)
{
  check_njacks_init();
  
  //open all files
  vector<string> file_list;
  for(auto &id : id_list) file_list.push_back(combine(template_path.c_str(),id));
  
  return read_conf_set_t(file_list,ntot_col,cols,nlines,verbosity);
}

djvec_t read_conf_set_t(const string &template_path,const range_t &range,size_t ntot_col,const vector<size_t> &cols,size_t nlines,bool verbosity)
{
  //get the list of existing files
  vector<size_t> id_list=get_existing_paths_in_range(template_path,range,verbosity);
  
  return read_conf_set_t(template_path,id_list,ntot_col,cols,nlines,verbosity);
}
