#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <boot.hpp>
#include <file.hpp>
#include <iostream>
#include <oper.hpp>

djvec_t read_conf_set_t(string template_path,size_t start,size_t each,size_t end,size_t ntot_col,vector<size_t> cols,int nlines)
{
  //basic checks
  check_njacks_init();
  if(end<start) CRASH("End=%d must be larger than Start=%d",end,start);
  if(each==0) CRASH("Each=0");
  
  //allocate files and open them
  size_t nconfs=(end+1-start)/each;
  vector<obs_file_t> files(nconfs,obs_file_t(ntot_col));
  for(size_t ind=0;ind<files.size();ind++)
    {
      files[ind].open(combine(template_path.c_str(),ind*each+start).c_str());
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
	  cout<<"size: "<<temp.size()<<endl;
	  cout<<raw_data[ind]<<temp<<endl;
	  
	  //if(temp.size()) raw_data[ind].insert(raw_data[ind].end(),temp.begin(),temp.end());
	}
      while(temp.size());
    }
  
  djvec_t out(raw_data);
  return out;
}
