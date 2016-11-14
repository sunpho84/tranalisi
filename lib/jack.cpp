#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_JACK
#include <jack.hpp>
#include <meas_vec.hpp>
#ifdef USE_OMP
 #include <omp.h>
#endif
#include <oper.hpp>
#include <tools.hpp>

void set_njacks(int ext_njacks)
{
  if(njacks==UNDEF_NJACKS) njacks=ext_njacks;
  else CRASH("Unbale to set njacks twice");
}

void check_njacks_init()
{if(njacks==UNDEF_NJACKS) CRASH("Set njacks before");}

djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,int nlines=1)
{
  //basic checks
  check_njacks_init();
  if(range.end<range.start) CRASH("End=%d must be larger than Start=%d",range.end,range.start);
  if(range.each==0) CRASH("Each=0");
  
  //allocate files and open them
  size_t nconfs=(range.end+1-range.start)/range.each;
  vector<obs_file_t> files(nconfs,obs_file_t(ntot_col));
  for(size_t ind=0;ind<files.size();ind++)
    {
      files[ind].open(combine(template_path.c_str(),ind*range.each+range.start).c_str());
      files[ind].set_col_view(cols);
    }
  
  //measure the length of the first file
  size_t length=files[0].length(nlines);
  cout<<"Total length: "<<length<<endl;
  
  //! output
  djvec_t data(length);
  
  //! raw data
  vector<vector<double>> raw_data(files.size());
#ifdef USE_OMP
 #pragma omp parallel for
#endif
  for(size_t ind=0;ind<files.size();ind++)
    {
#ifdef USE_OMP
      printf("Thread %d reading file %zu/%zu\n",omp_get_thread_num(),ind,files.size());
#else
      printf("Reading file %zu/%zu\n",ind,files.size());
#endif
      vector<double> temp;
      do
	{
	  temp=files[ind].read(nlines);
	  if(temp.size()) raw_data[ind].insert(raw_data[ind].end(),temp.begin(),temp.end());
	}
      while(temp.size());
    }
  cout<<"Finished reading"<<endl;
  
  return djvec_t(transpose(raw_data));
}
