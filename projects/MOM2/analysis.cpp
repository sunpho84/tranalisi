#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_ANALYSIS
 #include <MOM2/analysis.hpp>

#define ASSERT_COMPATIBLE_MEMBER(MEMBER)			\
  if(data(in1,ASSERT_PRESENT).MEMBER!=				\
     data(in1,ASSERT_PRESENT).MEMBER)				\
    CRASH("Impossible to combine, different member " #MEMBER)

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
  data(out,PRESENCE_NOT_NEEDED)=data(in1,ASSERT_PRESENT);
  data(out,PRESENCE_NOT_NEEDED).dir_path=out;
  
  for(auto &p : data(in2,ASSERT_PRESENT).get_all_tasks(data(out,PRESENCE_NOT_NEEDED)))
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

void sea_chir_extrap(const string out,const vector<string> &ens_list)
{
}

auto assert_ens_present(const string &key)
{
  auto f=_data.find(key);
  
  if(f==_data.end())
    {
      cout<<"Available keys: "<<endl;
      for(auto d : _data) cout<<" "<<d.first<<endl;
      CRASH("Unable to find the key %s",key.c_str());
    }
  
  return f;
}

void data_erase(const string &key)
{
  _data.erase(assert_ens_present(key));
}

perens_t& data(const string &key,const bool assert_present_flag)
{
  if(assert_present_flag)
    return assert_ens_present(key)->second;
  else
    return _data[key];
}

void list_ensembles()
{
  cout<<"Ensembles:"<<endl;
   for(auto &path : pars::ens)
     cout<<" "<<path<<endl;
}
