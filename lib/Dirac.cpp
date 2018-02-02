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

Dirac_t init_Gamma_onecol(const int *irow,const int *re,const int *im,size_t icol_row,size_t icol_col)
{
  //! matrix being built
  Dirac_t m(NSPINCOL,NSPINCOL);
  
  //! list of entries
  vector<Triplet<dcompl_t>> tr;
  
   for(size_t ispin=0;ispin<NSPIN;ispin++)
     tr.push_back(Triplet<dcompl_t>(ispin*NCOL+icol_col,irow[ispin]*NCOL+icol_row,dcompl_t(re[ispin],im[ispin])));
  
  m.setFromTriplets(tr.begin(),tr.end());
  
  return m;
}

namespace
{
  typedef int Gamma_data_t[3][4];
  
  Gamma_data_t Gamma_data[nGamma]=
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
  
  Gamma_data_t vit_Gamma_data[nGamma]=
    {{{0,1,2,3},{1,1,1,1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{1,1,-1,-1}},
     {{3,2,1,0},{1,-1,-1,1},{0,0,0,0}},
     {{2,3,0,1},{0,0,0,0},{1,-1,-1,1}},
     {{0,1,2,3},{1,1,-1,-1},{0,0,0,0}},
     {{2,3,0,1},{1,1,1,1},{0,0,0,0}},
     {{1,0,3,2},{0,0,0,0},{1,1,-1,-1}},
     {{1,0,3,2},{1,-1,-1,1},{0,0,0,0}},
     {{0,1,2,3},{0,0,0,0},{1,-1,-1,1}},
     {{2,3,0,1},{1,1,-1,-1},{0,0,0,0}},
     {{3,2,1,0},{0,0,0,0},{1,1,1,1}},
     {{3,2,1,0},{1,-1,1,-1},{0,0,-0,0}},
     {{2,3,0,1},{0,0,0,-0},{1,-1,1,-1}},
     {{1,0,3,2},{0,0,0,0},{1,1,1,1}},
     {{1,0,3,2},{1,-1,1,-1},{0,-0,0,0}},
     {{0,1,2,3},{0,-0,0,0},{1,-1,1,-1}}};
}

vector<Dirac_t> init_all_Gamma(const Gamma_data_t *data)
{
  //! returned list of Gamma
  vector<Dirac_t> out;
  
  for(size_t iGamma=0;iGamma<nGamma;iGamma++)
    out.push_back(init_Gamma(data[iGamma][0],data[iGamma][1],data[iGamma][2]));
  
  return out;
}

vector<Dirac_t> Gamma=init_all_Gamma(Gamma_data);

vector<Dirac_t> vit_Gamma=init_all_Gamma(vit_Gamma_data);

qprop_t convert_basis(const qprop_t &p,const Gamma_data_t *rec,const Gamma_data_t *pro)
{
  qprop_t pv;
  pv.Zero();
  for(size_t iG=0;iG<nGamma;iG++)
    for(size_t ic1=0;ic1<NCOL;ic1++)
      for(size_t ic2=0;ic2<NCOL;ic2++)
	pv+=(p*init_Gamma_onecol(pro[iG][0],pro[iG][1],pro[iG][2],ic1,ic2).adjoint()).trace()/4.0*
	  init_Gamma_onecol(rec[iG][0],rec[iG][1],rec[iG][2],ic1,ic2);
  
  return pv;
}

qprop_t convert_to_Vit_basis(const qprop_t &p)
{return convert_basis(p,vit_Gamma_data,Gamma_data);}

qprop_t convert_from_Vit_basis(const qprop_t &p)
{return convert_basis(p,Gamma_data,vit_Gamma_data);}

jqprop_t invert(const jqprop_t &in)
{
  jqprop_t out;
  for(size_t ijack=0;ijack<=njacks;ijack++) out[ijack]=in[ijack].inverse();
  return out;
}
