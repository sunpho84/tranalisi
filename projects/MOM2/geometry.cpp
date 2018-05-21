#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#ifdef USE_OMP
 #include <omp.h>
#endif

#include <fstream>
#include <iostream>

#define EXTERN_GEOMETRY
 #include <MOM2/geometry.hpp>

#include <MOM2/perens.hpp>

using namespace std;

double ph_mom[NDIM];

void perens_t::set_comp_list_of_moms(const string &mom_list_path,double filter_thresh)
{
  //slightly increment thresh to include border
  filter_thresh*=1+1e-10;
  
  //open the file to read momentum
  ifstream mom_file(mom_list_path);
  if(not mom_file.good()) CRASH("Unable to open %s",mom_list_path.c_str());
  
  //open the file to write filtered momentumm
  const string filt_path=dir_path+"/filt_mom.txt";
  ofstream filt_file(filt_path);
  if(not filt_file.good()) CRASH("Unable to open %s",filt_path.c_str());
  
  do
    {
      //temporary read the coords
      imom_t c;
      for(auto &ci : c) mom_file>>ci;
      //if coords good, store them
      if(mom_file.good())
	{
	  //store in any case
	  all_moms.push_back(c);
	  
	  //mark if filtered or not
	  const double discr=c.p(L).tilde().p4_fr_p22();
	  const bool filt=(discr<filter_thresh);
	  filt_moms.push_back(filt);
	  
	  filt_file<<c<<" = "<<c.p(L)<<" , discr: "<<discr<<" , filt: "<<filt<<endl;
	}
    }
  while(mom_file.good());
  
  //count the computed momenta
  ncomp_moms=all_moms.size();
  
  //count the number of momenta that passed the filter
  const size_t nthresh=count_if(filt_moms.begin(),filt_moms.end(),[](const bool &i){return i;});
  
  //print stats
  cout<<"Read "<<filt_moms.size()<<" momenta"<<endl;
  cout<<"NFiltered moms (p4/p2^2<"<<filter_thresh<<"): "<<nthresh<<endl;
}

void perens_t::write_comp_list_of_moms(const string &mom_list_path) const
{
  //open the file to write
  ofstream mom_file(mom_list_path);
  if(not mom_file.good()) CRASH("Unable to open %s",mom_list_path.c_str());
  
  for(const auto &l : linmoms)
    {
      const imom_t &c=all_moms[l[0]];
      for(auto &ci : c) mom_file<<ci<<" ";
      mom_file<<endl;
    }
}

size_t perens_t::get_mir_mom(size_t imom,size_t imir)
{
  coords_t cm=all_moms[imom];
  
  if(imir&1) cm[0]=-cm[0]-1;
  for(size_t mu=1;mu<NDIM;mu++) cm[mu]*=(1-2*((imir>>mu)&1));
  auto ret=find(all_moms.begin(),all_moms.end(),cm);
  if(ret==all_moms.end()) CRASH("searching imir=%zu of %zu",imom,imom);
  
  return distance(all_moms.begin(),ret);
}

void perens_t::set_ri_mom_moms()
{
  for(size_t imom=0;imom<ncomp_moms;imom++)
    if(filt_moms[imom])
      {
	const size_t linmom=linmoms.size();
	linmoms.push_back({imom});
	bilmoms.push_back({imom,linmom,linmom});
    }
}

void perens_t::set_smom_moms()
{
  //http://xxx.lanl.gov/pdf/0901.2599v2
  linmoms.clear();
  bilmoms.clear();
  
  const double tol=1e-10;
  for(size_t i=0;i<ncomp_moms;i++)
    if(filt_moms[i])
      {
	//get norm of p[i]
	p_t pi=all_moms[i].p(L);
	double pi2=pi.norm2();
	
	for(size_t j=0;j<ncomp_moms;j++)
	  if(filt_moms[j])
	    {
	      //get norm of p[j]
	      p_t pj=all_moms[j].p(L);
	      double pj2=pj.norm2();
	      
	      //check that norm of the two incoming vector is the same
	      if(2.0*fabs(pi2-pj2)<(pi2+pj2)*tol)
		{
		  //sum and get norm
		  imom_t momk;
		  for(size_t mu=0;mu<NDIM;mu++)
		    momk[mu]=all_moms[i][mu]-all_moms[j][mu];
		  double pk2=momk.p(L).norm2();
		  
		  //debug info
		  //cerr<<"pi2: "<<pi2<<" pk2: "<<pk2<<" "<<2.0*fabs(pi2-pk2)/(pi2+pk2)<<" "<<pi2<<" "<<pj2<<" "<<pk2<<endl;
		  
		  if(2.0*fabs(pi2-pk2)<(pi2+pk2)*tol)
		    {
		      //search in list
		      auto posk=find(all_moms.begin(),all_moms.end(),momk);
		      
		      //if not found, push into the list of glb_moms
		      if(posk==all_moms.end())
		      	{
		      	  posk=all_moms.end();
		      	  all_moms.push_back(momk);
		      	}
		      
		      constexpr bool debug=false;
		      
		       const size_t k=distance(all_moms.begin(),posk);
		      //inform
		      if(debug)
			cout<<"Found smom pair: "<<i<<all_moms[i]<<pi2<<" + "<<j<<all_moms[j]<<pj2<<" = "<<momk<<pk2<<endl;
		      vector<size_t> pos;
		      
		      //search in the linmoms: if found take the distance, otherwise add
		      for(const size_t ic : {i,j})
			{
			  if(debug)
			    cout<<"searching for "<<ic<<endl;
			  auto pos_ic=find(linmoms.begin(),linmoms.end(),array<size_t,1>{ic});
			  size_t d;
			  if(pos_ic==linmoms.end())
			    {
			      //the position will be the end
			      d=linmoms.size();
			      //include it
			      linmoms.push_back({ic});
			      if(debug)
				cout<<" not found"<<endl;
			    }
			  else
			    {
			      d=distance(linmoms.begin(),pos_ic);
			      if(debug)
				cout<<" found"<<endl;
			    }
			  
			  //add to the list
			  if(debug)
			    cout<<"Position: "<<d<<endl;
			  pos.push_back(d);
			}
		      
		      //store
		      bilmoms.push_back({k,pos[0],pos[1]});
		    }
		  // else
		  //    cout<<"Unable to find it"<<momk<<"="<<glb_moms[i]<<"+"<<glb_moms[j]<<endl;
		}
	    }
      }
  
  cout<<"Number of smom pairs: "<<bilmoms.size()<<endl;
}
