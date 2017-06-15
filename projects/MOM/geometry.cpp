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
  map<imom_t,vector<size_t>> equiv_imoms_map;
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
      equiv_imoms_map[cr].push_back(i);
    }
  
  //reserve the number of components different from 0 in each class
  Np_class.reserve(equiv_imoms.size());
  
  //trasform map into vector
  for(auto &mom_class : equiv_imoms_map)
    {
      imom_t repr=imoms[mom_class.second[0]]; //take the first real mom, as the key could not exist
      auto it=find(imoms.begin(),imoms.end(),repr);
      if(it==imoms.end()) CRASH("Something went wrong with %zu %zu %zu %zu",repr[0],repr[1],repr[2],repr[3]);
      equiv_imoms.push_back(make_pair(distance(imoms.begin(),it),mom_class.second));
      
      //compute the number of components different from 0 in each class
      size_t Np=0;
      for(size_t mu=0;mu<NDIM;mu++) Np+=(repr[mu]!=0 or fabs(ph_mom[mu])>1.0e-10);
      Np_class.push_back(Np);
    }
  
  //print stats
  cout<<"Found "<<equiv_imoms.size()<<" independent momenta "<<endl;
  
  //print details
  // auto Np=Np_class.begin();
  // for(auto &e : equiv_imoms)
  //   {
  //     for(auto &ei : imoms[e.first]) cout<<ei<<" ";
  //     cout<<"(Np="<<*(Np++)<<") ";
  //     cout<<"Equivalent to: "<<e.second.size()<<" moms:"<<endl;
  //     for(auto &eq : e.second)
  // 	{
  // 	  cout<<" ";
  // 	  for(auto &eqi : imoms[eq]) cout<<eqi<<" ";
  // 	  cout<<endl;
  // 	}
  //   }
}

void list_all_smom_pairs()
{
  size_t npairs=0;
  
  for(size_t i=0;i<imoms.size();i++)
    {
      p_t pi=imoms[i].p(L);
      double pi2=pi.norm2();
      for(size_t j=0;j<imoms.size();j++)
	{
	  p_t pj=imoms[j].p(L);
	  double pj2=pj.norm2();
	  
	  if(2.0*fabs(pi2-pj2)<(pi2+pj2)*1e-10)
	    {
	      p_t pk;
	      for(size_t mu=0;mu<NDIM;mu++) pk[mu]=pi[mu]+pj[mu];
	      double pk2=pk.norm2();
	      
	      cerr<<2.0*fabs(pi2-pk2)/(pi2+pk2)<<" "<<pi2<<" "<<pj2<<" "<<pk2<<endl;
	      
	      if(2.0*fabs(pi2-pk2)<(pi2+pk2)*1e-10)
		{
		  cout<<"Found smom pair: "<<i<<" "<<j<<endl;
		  npairs++;
		}
	    }
	}
    }
  
  cout<<"Number of smom pairs: "<<npairs<<endl;
}

vector<double> get_indep_pt2()
{
  vector<double> out;
  out.reserve(equiv_imoms.size());
  for(auto &mom_class : equiv_imoms) out.push_back(imoms[mom_class.first].p(L).tilde().norm2());
  
  return out;
}
