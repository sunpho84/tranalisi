#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <geometry.hpp>

//////////////// global variables ////////////////

//! read the input file
void parse_input(const string &name)
{
  raw_file_t input(name,"r");
  
  //! lattice size
  size_t Ls=input.read<size_t>("L");
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  L[0]=input.read<size_t>("T");
  
  //! list of momenta
  string mom_list_path=input.read<string>("MomList");
  get_list_of_moms(mom_list_path);
  get_class_of_equiv_moms();
  
  //print stats
  cout<<"Read "<<moms.size()<<" momenta"<<endl;
  cout<<"Found "<<equiv_moms.size()<<" independent momenta "<<endl;
  
  // //print details
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
  
  //! number of jacks
  const size_t ext_njacks=input.read<size_t>("NJacks");
  set_njacks(ext_njacks);
}

int main(int narg,char **arg)
{
  //read input file
  string name="input.txt";
  if(narg>=2) name=arg[1];
  parse_input(name);
  
  
  cout<<Gamma[5]<<endl;
  
  
  //a+=b;
  
  RowVectorXi adew;
  //no no no fai jack l'indice più interno, punto
  //cout<<a<<endl;
  cout<<Gamma[4]<<endl;
  
  //cout<<del<<endl;
  
  //testing complex jack
  jack_t<dcomplex> d,e;
  d+=e;

  
  return 0;
}
