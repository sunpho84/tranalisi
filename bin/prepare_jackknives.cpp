#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s [input]",arg[0]);
  
  //! initial time
  auto start=take_time();
  
  //! input file
  input_file_t input(arg[1]);
  
  //! the list of files
  size_t nfiles=input.read<size_t>("NFiles");
  vector<pair<string,string>> files(nfiles);
  for(auto &it : files)
    {
      it.first=input.read<string>("TemplateIn");
      it.second=input.read<string>("TemplateOut");
    }
  
  //! range of file
  range_t file_range;
  input.expect("ConfRange");
  file_range.start=input.read<size_t>();
  file_range.each=input.read<size_t>();
  file_range.end=input.read<size_t>();
  
  //! set the number of jackknives
  set_njacks(input.read<size_t>("NJacks"));
  
  //! time extension
  size_t T=input.read<size_t>("T");
  
  //! total number of columns
  size_t ntot_cols=input.read<size_t>("NTotCols");
  
  //! columns to read
  size_t ncols=input.read<size_t>("NCols");
  vector<size_t> cols(ncols);
  for(auto &it : cols) it=input.read<size_t>();
  
  //! list of filters to be made
  size_t nfilters=input.read<size_t>("NFilters");
  map<string,filter_t> filters;
  for(size_t ifilter=0;ifilter<nfilters;ifilter++)
    {
      string name=input.read<string>("Filter");
      size_t each=input.read<size_t>("each");
      size_t offset=input.read<size_t>("offset");
      filters[name]=filter_t(each*T*ncols,offset*T*ncols,T*ncols);
    }
  
  //loop on files
  for(auto &it : files)
    {
      djvec_t data=read_conf_set_t(it.first,file_range,ntot_cols,cols,T);
      
      //loop on filters
      for(auto &filter : filters)
	filter.second(data).bin_write(combine(it.second.c_str(),filter.first.c_str()));
    }

  cout<<"Total time: "<<elapsed_time(start)<<endl;
  
  return 0;
}
