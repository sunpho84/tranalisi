#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <grace.hpp>

grace_file_t &operator<<(grace_file_t &out,const vec_ave_err_t &data)
{
  out<<"@type xydy"<<endl;
  for(size_t i=0;i<data.size();i++) out<<i<<" "<<data[i]<<endl;
  out.new_set();
  return out;
}
