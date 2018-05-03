#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_PARS
 #include <MOM2/pars.hpp>

#include <MOM2/pr_bil.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/sigma.hpp>

void freeze_pars()
{
  pars::can_change_pars=false;
  
  lprop::set_ins();
  qprop::set_ins();
  jqprop::set_ins();
  
  sigma::set_ins();
  pr_bil::set_ins();
}
