#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <geometry.hpp>

//////////////// global variables ////////////////

//! list of momenta
vector<coords_t> moms;

//! read the input file
void parse_input(const string &name)
{
  raw_file_t input(name,"r");
  
  //! lattice size
  size_t Ls=input.read<size_t>("L");
  for(size_t mu=1;mu<NDIM;mu++) L[mu]=Ls;
  L[0]=input.read<size_t>("T");
  
  //! number of ranges
  const int nfft_ranges=input.read<int>("NFftRanges");
  
  //read all the ranges
  for(int irange=0;irange<nfft_ranges;irange++)
    {
      size_t L[2],T[2];
      L[0]=input.read<size_t>("L");
      L[1]=input.read<size_t>();
      T[0]=input.read<size_t>("T");
      T[1]=input.read<size_t>();
      
      //init the offset and width from range interval
      coords_t offs,width;
      offs[0]=T[0];
      width[0]=T[1]-T[0]+1;
      for(int i=1;i<NDIM;i++)
	{
	  offs[i]=L[0];
	  width[i]=L[1]-L[0]+1;
	}
      
      //put all momenta
      for(int vol=vol_of_lx(width),ifilt=0;ifilt<vol;ifilt++)
	{
	  //gets the coordinate in the filtering volume
	  coords_t c=offs+coord_of_lx(ifilt,width);
	  
	  //get mirrorized
	  for(int imir=0;imir<pow(2,NDIM);imir++)
	    {
	      coords_t cmir=get_mirrorized_site_coords(c,imir);
	      moms.push_back(cmir);
	    }
	}
    }
  
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
  
  jprop_t a,b;
  
  a.Zero();
  b=b.Identity();
  a+=b;
  
  a*=b;
  
  cout<<b<<endl;
  
  return 0;
}
