#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <ave_err.hpp>

ostream& operator<<(ostream &out,const ave_err_t &ae)
{
  if(ae.is_printable()) out<<ae.ave()<<" "<<ae.err();
  
  return out;
}
