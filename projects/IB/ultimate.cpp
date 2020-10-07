#include <common.hpp>

#define PROVIDE(A)				\
  const auto A=lat_par[input_an_id].A

int  main(int narg,char **arg)
{
  if(narg<2)
    CRASH("Error, use %s fileout",arg[0]);
  
  init_common_IB("/home/francesco/QCD/LAVORI/IB_MES_NF211/one_source/ultimate_input.txt");
  
  ofstream out(arg[1]);
  out.precision(16);
  out<<endl<<"# ml,B0,f0,ainv[0],ainv[1],ainv[2],Z[0],Z[1],Z[2]"<<endl<<endl;
  
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      out<<"# Branch: "<<input_an_id<<endl;
      
      PROVIDE(ml);
      PROVIDE(B0);
      PROVIDE(f0);
      PROVIDE(ainv);
      PROVIDE(Z);
      
      const vector<dboot_t> data{ml,B0,f0,ainv[0],ainv[1],ainv[2],Z[0],Z[1],Z[2]};
      
      out<<endl<<"# Averages"<<endl<<endl;
      for(auto& x : data)
	out<<x.ave()<<endl;
      
      out<<endl<<"# Cov Matrix"<<endl<<endl;
      
      for(auto& x : data)
	{
	  for(auto& y : data)
	    out<<cov(x,y)<<endl;
	  out<<endl;
	}
      out<<endl;
    }
  
  return 0;
}

