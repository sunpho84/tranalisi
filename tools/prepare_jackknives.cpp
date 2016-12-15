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
  map<string,pair<size_t,size_t>> filters;
  for(size_t ifilter=0;ifilter<nfilters;ifilter++)
    {
      string name=input.read<string>("Filter");
      size_t each=input.read<size_t>("each");
      size_t offset=input.read<size_t>("offset");
      filters[name]=make_pair(offset,each);
    }
  
  //loop on files
  for(auto &it : files)
    {
      auto start=take_time();
      djvec_t data=read_conf_set_t(it.first,file_range,ntot_cols,cols,T);
      cout<<"Time to read: "<<elapsed_time(start)<<endl;
      
      //loop on filters
      for(auto &filter : filters)
	{
	  size_t base_nel=T*ncols;
	  size_t offset=filter.second.first;
	  size_t each=filter.second.second;
	  size_t hw=data.size()/(base_nel*each);
	  
	  vec_filter(data,gslice(base_nel*offset,{hw,T*ncols},{each*T*ncols,1})).bin_write(combine(it.second.c_str(),filter.first.c_str()));
	  //((djvec_t)data[gslice(base_nel*offset,{hw,T*ncols},{each*T*ncols,1})]).bin_write(combine(it.second.c_str(),filter.first.c_str()));
	}
    }
  
  cout<<"Total time: "<<elapsed_time(start)<<endl;
  
  return 0;
}
