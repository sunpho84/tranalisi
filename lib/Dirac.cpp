#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_DIRAC
 #include <Dirac.hpp>

namespace
{
  const Gamma_data_t Gamma_data[nGamma]=
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
}

vector<Dirac_t> quaGamma=init_all_Gamma<NCOL>(Gamma_data);

vector<Dirac_t> lepGamma=init_all_Gamma<1>(Gamma_data);

jqprop_t invert(const jqprop_t &in)
{
  jqprop_t out;
  for(size_t ijack=0;ijack<=njacks;ijack++)
    out[ijack]=in[ijack].inverse();
  return out;
}
