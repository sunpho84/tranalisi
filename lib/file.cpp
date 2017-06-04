#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <iostream>

#include <file.hpp>

vector<size_t> get_existing_paths_in_range(const string &template_path,const range_t &range,bool verbosity)
{
  //basic checks
  if(range.end<range.start) CRASH("End=%d must be larger than Start=%d",range.end,range.start);
  if(range.each==0) CRASH("Each=0");
  
  //try to open all of them
  vector<size_t> id_list;
  for(size_t icheck_file=0;icheck_file<range.size();icheck_file++)
    {
      size_t id=range(icheck_file);
      string path=combine(template_path.c_str(),id);
      //cout<<"Considering file: "<<path<<endl;
      if(file_exists(path)) id_list.push_back(id);
      else if(verbosity) cout<<"Skipping unavailable file "<<path<<endl;
    }
  
  if(verbosity) cout<<"Found "<<id_list.size()<<" out of "<<range.size()<<endl;
  
  return id_list;
}
