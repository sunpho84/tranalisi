#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_COMMON
#include <common.hpp>

dboot_t read_boot(const raw_file_t &file)
{
  dboot_t out;
  for(size_t ib=0;ib<nboots;ib++) file.read(out[ib]);
  return out;
}

void init_common_IB(string ens_pars)
{
  set_njacks(15);
  
    raw_file_t file(ens_pars,"r");
  
  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].ainv[ibeta][nboots]=0;
      for(size_t iboot=0;iboot<nboots;iboot++)
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  {
	    file.read(lat_par[input_an_id].ainv[ibeta][iboot]);
	    lat_par[input_an_id].ainv[ibeta][nboots]+=lat_par[input_an_id].ainv[ibeta][iboot];
	  }
      for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].ainv[ibeta][nboots]/=nboots;
    }
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an/2;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[input_an_id].Z[ibeta][iboot]);
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=ninput_an/2;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(lat_par[input_an_id].Z[ibeta][iboot]);
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_total;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[input_an_id][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].B0=read_boot(file)/2.0;
}

//! perform the analysis according to eq.28
ave_err_t eq_28_analysis(const dbvec_t &v)
{
  ave_err_t ae;
  double sigma=0;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave();
      double e=v[i].err();
      ae.ave+=a;
      ae.err+=sqr(a);
      sigma+=sqr(e);
    }
  ae.ave/=v.size();
  ae.err/=v.size();
  sigma/=v.size();
  ae.err-=sqr(ae.ave);
  ae.err=sqrt(fabs(ae.err)+sigma);
  
  return ae;
}

