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

void set_glb_list_of_moms(const string &path,double thresh)
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
	  if(filt) glb_moms.push_back(c);
	  filt_moms.push_back(filt);
	}
    }
  while(mom_file.good());
  
  //print stats
  cout<<"Read "<<filt_moms.size()<<" momenta"<<endl;
  cout<<"NFiltered moms (p4/p2^2<"<<thresh<<"): "<<glb_moms.size()<<endl;
}

size_t get_mir_mom(size_t imom,size_t imir)
{
  coords_t cm=glb_moms[imom];
  
  if(imir&1) cm[0]=-cm[0]-1;
  for(size_t mu=1;mu<NDIM;mu++) cm[mu]*=(1-2*((imir>>mu)&1));
  auto ret=find(glb_moms.begin(),glb_moms.end(),cm);
  if(ret==glb_moms.end()) CRASH("searching imir=%zu of %zu",imom,imom);
  
  return distance(glb_moms.begin(),ret);
}
// void list_all_smom_pairs()
// {
//   size_t npairs=0;
  
//   for(size_t i=0;i<imoms.size();i++)
//     {
//       p_t pi=imoms[i].p(L);
//       double pi2=pi.norm2();
//       for(size_t j=0;j<imoms.size();j++)
// 	{
// 	  p_t pj=imoms[j].p(L);
// 	  double pj2=pj.norm2();
	  
// 	  if(2.0*fabs(pi2-pj2)<(pi2+pj2)*1e-10)
// 	    {
// 	      p_t pk;
// 	      for(size_t mu=0;mu<NDIM;mu++) pk[mu]=pi[mu]+pj[mu];
// 	      double pk2=pk.norm2();
	      
// 	      cerr<<2.0*fabs(pi2-pk2)/(pi2+pk2)<<" "<<pi2<<" "<<pj2<<" "<<pk2<<endl;
	      
// 	      if(2.0*fabs(pi2-pk2)<(pi2+pk2)*1e-10)
// 		{
// 		  cout<<"Found smom pair: "<<i<<" "<<j<<endl;
// 		  npairs++;
// 		}
// 	    }
// 	}
//     }
  
//   cout<<"Number of smom pairs: "<<npairs<<endl;
// }

// vector<double> get_pt2()
// {
//   vector<double> out;
//   out.reserve(imoms.size());
//   for(auto &mom : imoms) out.push_back(mom.p(L).tilde().norm2());
  
//   return out;
// }

// vector<double> get_indep_pt2()
// {
//   vector<double> out;
//   out.reserve(equiv_imoms.size());
//   for(auto &mom_class : equiv_imoms) out.push_back(imoms[mom_class.first].p(L).tilde().norm2());
  
//   return out;
// }

// djvec_t average_equiv_moms(const djvec_t &in,const index_t &oth_ind_mom_ind,const index_t &oth_mom_ind)
// {
//   djvec_t out(oth_ind_mom_ind.max());
// #pragma omp parallel for
//   for(size_t ioth_ind_mom=0;ioth_ind_mom<oth_ind_mom_ind.max();ioth_ind_mom++)
//     {
//       const vector<size_t> oth_ind_mom_comps=oth_ind_mom_ind(ioth_ind_mom);
//       const size_t ind_mom=oth_ind_mom_comps.back();
      
//       //reset
//       out[ioth_ind_mom]=0.0;
      
//       //loop on equivalent moms
//       auto &imom_class=equiv_imoms[ind_mom];
//       for(size_t imom : imom_class.second)
// 	{
// 	  //copy all components
// 	  vector<size_t> oth_mom_comps=oth_ind_mom_comps;
// 	  //set the last component to imom
// 	  oth_mom_comps.back()=imom;
// 	  //take index
// 	  size_t ioth_mom=oth_mom_ind(oth_mom_comps);
// 	  out[ioth_ind_mom]+=in[ioth_mom];
// 	}
      
//       //normalize
//       out[ioth_ind_mom]/=imom_class.second.size();
//     }
  
//   return out;
// }
