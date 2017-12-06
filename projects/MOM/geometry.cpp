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
