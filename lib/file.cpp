#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <file.hpp>
#include <iostream>
#include <oper.hpp>

djvec_t read_conf_set_t(string template_path,range_t range,size_t ntot_col,vector<size_t> cols,int nlines)
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
  
  //! raw data
  vector<vector<double>> raw_data(files.size());
  for(size_t ind=0;ind<files.size();ind++)
    {
      vector<double> temp;
      do
	{
	  //read
	  temp=files[ind].read(nlines);
	  //cout<<raw_data[ind]<<temp<<endl;
	  
	  //append
	  if(temp.size()) raw_data[ind].insert(raw_data[ind].end(),temp.begin(),temp.end());
	}
      while(temp.size());
    }
  
  return djvec_t(transpose(raw_data));;
}
