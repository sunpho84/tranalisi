#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <Dirac.hpp>

using namespace std;

#define PROVIDE_DECRYPTER						\
  									\
  inline auto decrypt(const string &str) -> decltype(decr.find(str)->second) \
  {									\
    auto act_key=decr.find(str);					\
    if(act_key==decr.end())						\
      {									\
	cout<<"Available keys: "<<endl;					\
	for(const auto &p : decr)					\
	  cout<<" "<<p.first<<endl;					\
	CRASH("Unable to decrypt %s",str.c_str());			\
      }									\
    return act_key->second;						\
  }									\
  /* repeat to swallow semicolon*/					\
  void decrypt()

#endif
