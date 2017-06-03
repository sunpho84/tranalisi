#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_GEOMETRY
#include <geometry.hpp>

#include <tools.hpp>
#include <fstream>

using namespace std;

void get_list_of_moms(const string &path)
{
  //open the file
  ifstream mom_file(path);
  if(not mom_file.good()) CRASH("Unable to open %s",path.c_str());
  do
    {
      //temporary read the coords
      mom_t c;
      for(auto &ci : c) mom_file>>ci;
      //if coords good, store them
      if(mom_file.good()) moms.push_back(c);
    }
  while(mom_file.good());
}

void get_class_of_equiv_moms()
{
  for(size_t i=0;i<moms.size();i++)
    {
      //get representative
      mom_t cr;
      for(size_t mu=0;mu<NDIM;mu++) cr[mu]=abs(moms[i][mu]);
      sort(&cr[1],cr.end());
      //store the index to equvalents
      equiv_moms[cr].push_back(i);
    }
}
