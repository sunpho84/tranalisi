#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_PARS
 #include <MOM2/pars.hpp>

#include <MOM2/sigma.hpp>

void freeze_pars()
{
  pars::can_change_pars=false;
  
  set_sigma_ins();
}
