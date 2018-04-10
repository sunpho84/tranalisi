#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

const size_t T=32;
const size_t nr=2;
const size_t nmel=5;
enum mel_t{V0P5,V1P5,V2P5,V3P5,P5P5};
const size_t nm=3;
double m_list[]={0.10,0.20,0.30};
const size_t nth=5;
double th_list[]={0.00,-0.10,0.10,-0.20,0.20};
const string base_run="/marconi_work/INF17_lqcd123_0/sanfo/STEP_SCALING/L16_T32_beta6.61/Semileptonic/";

djvec_t retrive_data(const string &binfile,const string &template_path,const index_t &ind)
{
  const range_t range{900,100,40800};
  const size_t ntot_col=2;
  const vector<size_t> cols={0,1};
  
  cout<<"npts: "<<ind.max()<<endl;
  
  djvec_t data;
  if(file_exists(binfile))
    {
      cout<<"Binary file \""<<binfile<<"\" exists, loading"<<endl;
      
      data.resize(ind.max());
      data.bin_read(binfile);
    }
  else
    {
      cout<<"Binary file \""<<binfile<<"\" does not exist, reading from txt confs, \""<<template_path<<"\""<<endl;
      
      data=read_conf_set_t(template_path,range,ntot_col,cols);
      if(data.size()!=ind.max()) CRASH("loaded binary file \"%c\" data has size %zu, expected %zu",binfile.c_str(),data.size(),ind.max());
      data.bin_write(binfile);
    }
  
  return data;
}

djvec_t get_2pts(const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw,const mel_t imel)
{
  static const index_t i2pts_ind({{"mbw",nm},{"rbw",nr},{"mfw",nm},{"rfw",nr},{"thbw",nth},{"thfw",nth},{"mel",nmel},{"T",T},{"ri",2}});
  static djvec_t data_2pts=retrive_data("2pts.dat",base_run+"out/%05d/mes_contr_2pts",i2pts_ind);
  
  int par;
  size_t ri;
  
  switch(imel)
    {
    case P5P5:
      par=+1;
      ri=0;
      break;
    default:
      par=0;
      ri=0;
      CRASH("Why!");
    }
  
  djvec_t out(T);
  
  for(size_t t=0;t<T;t++)
    for(size_t r=0;r<nr;r++)
      {
	size_t i=i2pts_ind({imbw,r,imfw,r,ithbw,ithfw,imel,t,ri});
	out[t]+=data_2pts[i];
      }
  
  return out.symmetrized(par)/nr;
}

djvec_t get_2pts_P5P5(const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw)
{
  return get_2pts(imbw,imfw,ithbw,ithfw,P5P5);
}

djvec_t get_3pts(const size_t imspec,const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw,const mel_t imel)
{
  static const index_t i3pts_ind({{"ms",nm},{"rs",nr},{"m1",nm},{"r1",nr},{"m0",nm},{"r0",nr},{"th1",nth},{"th0",nth},{"el",nmel},{"T",T},{"ri",2}});
  static djvec_t data_3pts=retrive_data("3pts.dat",base_run+"out/%05d/mes_contr_3pts",i3pts_ind);
  
  int par;
  size_t ri;
  
  switch(imel)
    {
    case V0P5:
      par=-1;
      ri=0;
      break;
    case V1P5:
    case V2P5:
    case V3P5:
      par=+1;
      ri=1;
      break;
    default:
      par=0;
      ri=0;
      CRASH("Why!");
    }
  
  djvec_t out(T);
  
  for(size_t t=0;t<T;t++)
    for(size_t r=0;r<nr;r++)
      {
	size_t i=i3pts_ind({imspec,r,imbw,!r,imfw,r,ithbw,ithfw,imel,t,ri});
	out[t]+=data_3pts[i];
      }
  
  return out.symmetrized(par)/nr;
}

djvec_t get_3pts_V0P5(const size_t imspec,const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw)
{
  return get_3pts(imspec,imbw,imfw,ithbw,ithfw,V0P5);
}

djvec_t get_3pts_VKP5(const size_t imspec,const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw)
{
  return (get_3pts(imspec,imbw,imfw,ithbw,ithfw,V1P5)+
	  get_3pts(imspec,imbw,imfw,ithbw,ithfw,V2P5)+
	  get_3pts(imspec,imbw,imfw,ithbw,ithfw,V3P5))
    /3.0;
}

djvec_t get_3pts_VKP5_impr(const size_t imspec,const size_t imbw,const size_t imfw,const size_t ithbw,const size_t ithfw)
{
  int th_opp[]={0,2,1,4,3};
  
  return (get_3pts_VKP5_impr(imspec,imbw,imfw,ithbw,ithfw)+
	  get_3pts_VKP5_impr(imspec,imbw,imfw,th_opp[ithbw],th_opp[ithfw]))/2.0;
}

int main()
{
  set_njacks(15);
  
  cout<<get_2pts_P5P5(0,0,0,0).ave_err()<<endl;
  cout<<get_3pts_V0P5(0,0,0,0,0).ave_err()<<endl;
  cout<<get_3pts_VKP5(0,0,0,1,2).ave_err()<<endl;
  cout<<get_3pts_VKP5(0,0,0,2,1).ave_err()<<endl;
  cout<<get_3pts_VKP5_impr(0,0,0,2,1).ave_err()<<endl;
  
  return 0;
}
