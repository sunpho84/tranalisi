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

void set_list_of_moms(const string &path,double thresh)
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
      if(mom_file.good())
	{
	  bool filt=(c.p(L).tilde().p4_fr_p22()<thresh);
	  if(filt) imoms.push_back(c);
	  filt_moms.push_back(filt);
	}
    }
  while(mom_file.good());
  
  //print stats
  cout<<"Read "<<filt_moms.size()<<" momenta"<<endl;
  cout<<"NFiltered moms (p4/p2^2<"<<thresh<<"): "<<imoms.size()<<endl;
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
      double p0=fabs(ph_mom[0]);
      if(fabs(p0)>1e-10 and fabs(p0-0.5)>1e-10) CRASH("phase on momentum 0 cannot be %lg",ph_mom[0]);
      if(fabs(p0-0.5)<=1e-10 and cr[0]<0) cr[0]=-cr[0]-1;
      
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
  ofstream out("equiv_moms.txt");
  out<<"Found "<<equiv_imoms.size()<<" independent momenta "<<endl;
  
  //print details
  for(auto &e : equiv_imoms)
    {
      out<<e.first<<" , p2hat: "<<imoms[e.first].p(L).tilde().norm2()<<endl;
      out<<" Equivalent to: "<<e.second.size()<<" moms: "<<endl;
      for(auto &eq : e.second)
  	{
  	  out<<"  "<<eq;
	  //components
  	  out<<"={";
	  for(size_t mu=0;mu<NDIM;mu++)
	    {
	      if(mu) out<<",";
	      out<<imoms[eq][mu];
	    }
	  out<<"}";
	  //phat
  	  out<<"={";
	  for(size_t mu=0;mu<NDIM;mu++)
	    {
	      if(mu) out<<",";
	      out<<imoms[eq].p(L).hat()[mu];
	    }
	  out<<"}"<<endl;
  	}
    }
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

djvec_t average_equiv_moms(const djvec_t &in,const index_t &oth_ind_mom_ind,const index_t &oth_mom_ind)
{
  djvec_t out(oth_ind_mom_ind.max());
#pragma omp parallel for
  for(size_t ioth_ind_mom=0;ioth_ind_mom<oth_ind_mom_ind.max();ioth_ind_mom++)
    {
      const vector<size_t> oth_ind_mom_comps=oth_ind_mom_ind(ioth_ind_mom);
      const size_t ind_mom=oth_ind_mom_comps.back();
      
      //reset
      out[ioth_ind_mom]=0.0;
      
      //loop on equivalent moms
      auto &imom_class=equiv_imoms[ind_mom];
      for(size_t imom : imom_class.second)
	{
	  //copy all components
	  vector<size_t> oth_mom_comps=oth_ind_mom_comps;
	  //set the last component to imom
	  oth_mom_comps.back()=imom;
	  //take index
	  size_t ioth_mom=oth_mom_ind(oth_mom_comps);
	  out[ioth_ind_mom]+=in[ioth_mom];
	}
      
      //normalize
      out[ioth_ind_mom]/=imom_class.second.size();
    }
  
  return out;
}
