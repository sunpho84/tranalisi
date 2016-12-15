#include <tranalisi.hpp>

using namespace std;

int main(int narg,char **arg)
{
  if(narg<4) CRASH("Use %s out c1 file1 [c2 file2 ... ]",arg[0]);
  if(narg%2) CRASH("Specify one coefficient for each file");

  //read the number of doubles and initialize the output vector
  size_t ndoubles=file_size(arg[2])/sizeof(double);
  vector<double> out(ndoubles,0.0);
  
  //loop over in files
  size_t nfiles_in=(narg-2)/2;
  for(size_t ifile_in=0;ifile_in<nfiles_in;ifile_in++)
    {
      //convert the coefficient
      double c;
      if(sscanf(arg[2+ifile_in*2],"%lg",&c)!=1) CRASH("Unable to convert to double \"%s\"",arg[2+ifile_in*2]);
      
      //open the file and read it
      raw_file_t files_in(arg[3+ifile_in*2],"r");
      for(size_t idouble=0;idouble<ndoubles;idouble++)
	{
	  double d;
	  files_in.bin_read(d);
	  out[idouble]+=c*d;
	}
    }
  
  //open the output file
  raw_file_t file_out(arg[1],"w");
  for(size_t idouble=0;idouble<ndoubles;idouble++) file_out.bin_write(out[idouble]);
  
  return 0;
}
