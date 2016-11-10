#include <tranalisi.hpp>

int main(int narg,char **arg)
{
  input_file_t input("/tmp/input");
  
  //read the list of files
  vector<pair<string,string>> file_list;
  int nfiles=input.read<int>("NFiles");
  for(int ifile=0;ifile<nfiles;ifile++)
    {
      
      string templ_in=input.read<string>("TemplateIn");
      string templ_out=input.read<string>("TemplateOut");
      //cout<<input.read<string>("TemplateOut")<<endl;
      file_list.push_back(make_pair(templ_in,templ_out));
    }
  
  return 0;
}
