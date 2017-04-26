#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif
#include <meas_vec.hpp>

djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,size_t nlines)
{
  //basic checks
  check_njacks_init();
  if(range.end<range.start) CRASH("End=%d must be larger than Start=%d",range.end,range.start);
  if(range.each==0) CRASH("Each=0");
  
  //compute nfiles
  size_t navail_files=(range.end-range.start)/range.each+1;
  
  //try to open all of them
  vector<obs_file_t> files;
  for(size_t icheck_file=0;icheck_file<navail_files;icheck_file++)
    {
      string path=combine(template_path.c_str(),icheck_file*range.each+range.start);
      //cout<<"Considering file: "<<path<<endl;
      if(file_exists(path)) files.push_back(obs_file_t(path,ntot_col,cols));
      else cout<<"Skipping unavailable file "<<path<<endl;
    }
  
  //trim
  cout<<"Opened "<<files.size()<<" out of "<<navail_files<<endl;
  size_t clust_size=files.size()/njacks;
  size_t nfiles=clust_size*njacks;
  files.resize(nfiles);
  cout<<"Trimmed to "<<nfiles<<", clust_size="<<clust_size<<endl;
  
  //measure the length of the first file
  size_t length=files[0].length(nlines);
  cout<<"Total length (in multiple of nlines="<<nlines<<"): "<<length<<endl;
  
  //check that the length is a multiple of ncols*nlines
  size_t block_nentr=nlines*cols.size();
  size_t nblocks=length/block_nentr;
  if(length!=nblocks*block_nentr) CRASH("Total length is not a multiple of ncols*nlines=%zu",block_nentr);
  
  //! output
  djvec_t data(length);
  
#pragma omp parallel for
  for(size_t ijack=0;ijack<njacks;ijack++)
    for(size_t ifile=ijack*clust_size;ifile<(ijack+1)*clust_size;ifile++)
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
  cout<<"Finished reading"<<endl;
  
  //clusterize each entry
  {
    auto start=take_time();
    for(auto &d : data) d.clusterize(clust_size);
    cout<<elapsed_time(start)<<" to clusterize"<<endl;
  }
  
  return data;
}
