#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <boot.hpp>
#include <oper.hpp>

ostream& operator<<(ostream &out,const ave_err_t &ae)
{
  if(!isnan(ae.ave) && !isnan(ae.err)) out<<ae.ave<<" "<<ae.err;
  
  return out;
}
