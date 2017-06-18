#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#define EXTERN_GEOMETRY
 #include <geometry.hpp>

#include <tools.hpp>
#include <fstream>
#include <iostream>

using namespace std;

double ph_mom[NDIM]={0.5,0,0,0};

void set_list_of_moms(const string &path)
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

size_t get_mir_mom(size_t imom,size_t imir)
{
  coords_t cm=imoms[imom];
  
  if(imir&1) cm[0]=-cm[0]-1;
  for(size_t mu=1;mu<NDIM;mu++) cm[mu]*=(1-2*((imir>>mu)&1));
  auto ret=find(imoms.begin(),imoms.end(),cm);
  if(ret==imoms.end()) CRASH("searching imir=%zu of %zu",imom,imom);
  
  return distance(imoms.begin(),ret);
}

void set_class_of_equiv_moms()
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
  
  //trasform map into vector
  for(auto &mom_class : equiv_imoms_map)
    {
      imom_t repr=imoms[mom_class.second[0]]; //take the first real mom, as the key could not exist
      auto it=find(imoms.begin(),imoms.end(),repr);
      if(it==imoms.end()) CRASH("Something went wrong with %zu %zu %zu %zu",repr[0],repr[1],repr[2],repr[3]);
      equiv_imoms.push_back(make_pair(distance(imoms.begin(),it),mom_class.second));
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

vector<double> get_pt2()
{
  vector<double> out;
  out.reserve(imoms.size());
  for(auto &mom : imoms) out.push_back(mom.p(L).tilde().norm2());
  
  return out;
}

vector<double> get_indep_pt2()
{
  vector<double> out;
  out.reserve(equiv_imoms.size());
  for(auto &mom_class : equiv_imoms) out.push_back(imoms[mom_class.first].p(L).tilde().norm2());
  
  return out;
}

vector<double> get_filtered_pt2()
{
  vector<double> out;
  out.reserve(iequiv_mom_of_ifilt.size());
  for(auto &ieq : iequiv_mom_of_ifilt)
    out.push_back(imoms[equiv_imoms[ieq].first].p(L).tilde().norm2());
  
  return out;
}

djvec_t average_equiv_moms(const djvec_t &in)
{
  djvec_t out(equiv_imoms.size());
#pragma omp parallel for
  for(size_t ind_mom=0;ind_mom<equiv_imoms.size();ind_mom++)
    {
      //reset
      out[ind_mom]=0.0;
      
      //loop on equivalent moms
      auto &imom_class=equiv_imoms[ind_mom];
      for(size_t imom : imom_class.second) out[ind_mom]+=in[imom];
      
      //normalize
      out[ind_mom]/=imom_class.second.size();
    }
  
  return out;
}

void set_filtered_moms(const double thresh)
{
  for(size_t ieq=0;ieq<equiv_imoms.size();ieq++)
    if(imoms[equiv_imoms[ieq].first].p(L).tilde().p4_fr_p22()<thresh)
      iequiv_mom_of_ifilt.push_back(ieq);
  cout<<"NFiltered moms (p4/p2^2<"<<thresh<<"): "<<iequiv_mom_of_ifilt.size()<<endl;
}

djvec_t get_filtered_moms(const djvec_t &in)
{
  djvec_t out(iequiv_mom_of_ifilt.size());
  for(size_t ifilt=0;ifilt<iequiv_mom_of_ifilt.size();ifilt++)
    out[ifilt]=in[iequiv_mom_of_ifilt[ifilt]];
  
  return out;
}
