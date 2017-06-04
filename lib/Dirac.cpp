#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_DIRAC
#include <Dirac.hpp>

Dirac_t init_Gamma(const int *irow,const int *re,const int *im)
{
  //! matrix being built
  Dirac_t m(NSPINCOL,NSPINCOL);
  
  //! list of entries
  vector<Triplet<dcompl_t>> tr;
  
  for(size_t ispin=0;ispin<NSPIN;ispin++)
    for(size_t icol=0;icol<NCOL;icol++)
      tr.push_back(Triplet<dcompl_t>(ispin*NCOL+icol,irow[ispin]*NCOL+icol,dcompl_t(re[ispin],im[ispin])));
  
  m.setFromTriplets(tr.begin(),tr.end());
  
  return m;
}

vector<Dirac_t> init_all_Gamma()
{
  //! returned list of Gamma
  vector<Dirac_t> out;
  
  int data[nGamma][3][4]=
    {{{0,1,2,3},{1,1,1,1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{-1,-1,1,1}},
     {{3,2,1,0},{-1,1,1,-1},{0,0,0,0}},
     {{2,3,0,1},{0,0,0,0},{-1,1,1,-1}},
     {{2,3,0,1},{-1,-1,-1,-1},{0,0,0,0}},
     {{0,1,2,3},{1,1,-1,-1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{1,1,1,1}},
     {{3,2,1,0},{1,-1,1,-1},{0,0,0,0}},
     {{2,3,0,1},{0,0,0,0},{1,-1,1,-1}},
     {{2,3,0,1},{1,1,-1,-1},{0,0,0,0}},
     {{1,0,3,2},{0,0,0,0},{-1,-1,1,1}},
     {{1,0,3,2},{-1,1,1,-1},{0,0,0,0}},
     {{0,1,2,3},{0,0,0,0},{-1,1,1,-1}},
     {{1,0,3,2},{0,0,0,0},{1,1,1,1}},
     {{1,0,3,2},{1,-1,1,-1},{0,0,0,0}},
     {{0,1,2,3},{0,0,0,0},{1,-1,1,-1}}};
  
  for(size_t iGamma=0;iGamma<nGamma;iGamma++)
    out.push_back(init_Gamma(data[iGamma][0],data[iGamma][1],data[iGamma][2]));
  
  return out;
}

vector<Dirac_t> Gamma=init_all_Gamma();
