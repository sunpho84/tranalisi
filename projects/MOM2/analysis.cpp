#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_ANALYSIS
 #include <MOM2/analysis.hpp>

#define ASSERT_COMPATIBLE_MEMBER(MEMBER)\
  if(data(in1).MEMBER!=data(in1).MEMBER)			\
    CRASH("Impossible to average, different member " #MEMBER)

void assert_compatible(const string in1,const string in2)
{
  for(size_t mu=0;mu<NDIM;mu++) ASSERT_COMPATIBLE_MEMBER(L[mu]);
  ASSERT_COMPATIBLE_MEMBER(nm);
  ASSERT_COMPATIBLE_MEMBER(nr);
  ASSERT_COMPATIBLE_MEMBER(linmoms);
  ASSERT_COMPATIBLE_MEMBER(bilmoms);
}

void average(const string out,const string in1,const string in2)
{
  assert_compatible(in1,in2);
  
  pars::ens.push_back(out);
  data(out)=data(in1);
  
  for(auto &p : data(in2).get_all_tasks(data(out)))
    {
      djvec_t &out=*p.out;
      const djvec_t &in=*p.in;
      
      out+=in;
      out/=2.0;
    }
  
  //remove from the list
  for(auto in : {in1,in2})
    {
      data_erase(in);
      pars::ens.erase(find(pars::ens.begin(),pars::ens.end(),in));
    }
}

void data_erase(const string &key)
{
  auto f=_data.find(key);
  
  if(f==_data.end())
    {
      cout<<"Available keys: "<<endl;
      for(auto d : _data) cout<<" "<<d.first<<endl;
      CRASH("Unable to find the key %s",key.c_str());
    }
  
  _data.erase(f);
}

perens_t& data(const string &key)
{
  return _data[key];
}

void list_ensembles()
{
  cout<<"Ensembles:"<<endl;
   for(auto &path : pars::ens)
     cout<<" "<<path<<endl;
}
