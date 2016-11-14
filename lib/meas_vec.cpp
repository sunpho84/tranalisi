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
  
  //compute nfiles and clust size
  size_t navail_files=(range.end+1-range.start)/range.each;
  size_t clust_size=navail_files/njacks;
  size_t nfiles=clust_size*njacks;
  if(navail_files!=nfiles) cout<<"Reducing nfiles from available "<<navail_files<<" to "<<nfiles<<endl;
  
  //allocate files and open them
  vector<obs_file_t> files(nfiles,obs_file_t(ntot_col));
  size_t icheck_file=0,ifile=0;
  do
    {
      //! full path
      string path=combine(template_path.c_str(),icheck_file*range.each+range.start);
      
      //if file exists open it, otherwise skip it
      if(file_exists(path))
	{
	  files[ifile].open(path);
	  files[ifile].set_col_view(cols);
	  ifile++;
	}
      else cout<<"Skipping unavailable file "<<path<<endl;
      icheck_file++;
    }
  while(icheck_file<navail_files && ifile<nfiles);
  
  //check that all files have been opened
  if(ifile<nfiles) CRASH("Unable to open all files, opened %zu instead of %zu",ifile,nfiles);
  
  //measure the length of the first file
  size_t length=files[0].length(nlines);
  cout<<"Total length: "<<length<<endl;
  
  //check that the length is a multiple of ncols*nlines
  size_t block_nentr=nlines*cols.size();
  size_t nblocks=length/block_nentr;
  if(length!=nblocks*block_nentr) CRASH("Total length is not a multiple of ncols*nlines=%zu",block_nentr);
  
  //! output
  djvec_t data(length);
  
#pragma omp parallel for
  for(size_t ijack=0;ijack<njacks;ijack++)
    for(size_t ifile=ijack*clust_size;ifile<=(ijack+1)*clust_size;ifile++)
      {
#ifdef USE_OMP
	cout<<"Thread "<<omp_get_thread_num()<<"/"<<omp_get_num_threads();
#endif
	cout<<"Reading file "<<ifile<<"/"<<files.size()<<endl;
	
	//read all blocks
	for(size_t iblock=0;iblock<nblocks;iblock++)
	  {
	    vector<double> temp=files[icheck_file].read(nlines);
	    if(temp.size()!=block_nentr) CRASH("Error reading file %zu, iblock %zu",ifile,iblock);
	    
	    //copy
	    for(size_t ientr=0;ientr<block_nentr;ientr++) data[ientr+block_nentr*iblock][ijack]+=temp[ientr];
	  }
    }
  cout<<"Finished reading"<<endl;
  
  //clusterize each entry
  {
    START_TIME();
    for(auto &d : data) d.clusterize(clust_size);
    cout<<ELAPSED_TIME()<<" to clusterize"<<endl;
  }
  
  return data;
}
