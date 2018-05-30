#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_SCAN
 #include <scan.hpp>

namespace scan
{
  void print_input_files_stack()
  {
    cout<<"Current stack of input files:"<<endl;
    for(auto& p : fin) cout<<"   "<<p.second<<endl;
  }
  
  void include(const string &path)
  {
    cout<<"Opening \""<<path<<"\""<<endl;
    
    //open
    ifstream file(path);
    
    //check failure
    if(not file.good())
      {
	print_input_files_stack();
	CRASH("Failed to open %s",path.c_str());
      }
    fin.push_back(make_pair(move(file),path));
  }
  
  void scan(const string &path)
  {
    include(path);
    
    init_scanner();
    parser_parse(nullptr);
    destroy_scanner();
  }
  
  void lex(char* buf,int &result,int max_size)
  {
    result=0;
    
    do
      {
	//point to last opened file
	ifstream& f=fin.back().first;
	
	//read
	f.get(*buf);
	
	//close file if end reached
	if(f.eof())
	  {
	    print_input_files_stack();
	    fin.pop_back();
	  }
	//check what read otherwise
	else
	  result=f.good();
      }
    //exit loop when no more file can be read or result is valid
    while(fin.size()>=1 and result==0);
  }
}
