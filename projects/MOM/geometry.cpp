#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_GEOMETRY
#include <geometry.hpp>

#include <tools.hpp>
#include <fstream>
#include <iostream>

using namespace std;

double ph_mom[NDIM]={0.5,0,0,0};

void get_list_of_moms(const string &path)
{
  //open the file
  ifstream mom_file(path);
  if(not mom_file.good()) CRASH("Unable to open %s",path.c_str());
  do
    {
      //temporary read the coords
      imom_t c;
      for(auto &ci : c) mom_file>>ci;
      //if coords good, store them
      if(mom_file.good()) imoms.push_back(c);
    }
  while(mom_file.good());
  
  //print stats
  cout<<"Read "<<imoms.size()<<" momenta"<<endl;
}

void get_class_of_equiv_moms()
{
  for(size_t i=0;i<imoms.size();i++)
    {
      //get representative
      imom_t cr;
      
      //decide time component
      cr[0]=imoms[i][0];
      if(cr[0]<0) cr[0]=-cr[0]-1;
      
      //decide space componnents
      for(size_t mu=1;mu<NDIM;mu++) cr[mu]=abs(imoms[i][mu]);
      sort(&cr[1],cr.end());
      
      //store the index to equvalents
      equiv_imoms[cr].push_back(i);
    }
  
  //print stats
  cout<<"Found "<<equiv_imoms.size()<<" independent momenta "<<endl;
  
  //print details
  // for(auto &e : equiv_moms)
  //   {
  //     for(auto &ei : e.first) cout<<ei<<" ";
  //     cout<<"Equivalent to: "<<e.second.size()<<" moms:"<<endl;
  //     for(auto &eq : e.second)
  // 	{
  // 	  cout<<" ";
  // 	  for(auto &eqi : moms[eq]) cout<<eqi<<" ";
  // 	  cout<<endl;
  // 	}
  //   }
}

vector<double> get_indep_pt2()
{
  vector<double> out;
  out.reserve(equiv_imoms.size());
  for(auto &mom_class : equiv_imoms) out.push_back(mom_class.first.p(L).tilde().norm2());
  
  return out;
}
