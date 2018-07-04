#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_PARS
 #include <MOM2/pars.hpp>

#include <MOM2/pr_bil.hpp>
#include <MOM2/pr_meslep.hpp>
#include <MOM2/prop.hpp>
#include <MOM2/sigma.hpp>

void freeze_pars()
{
  if(pars::can_change_pars)
    {
      if(pars::use_QED!=1 and (pars::use_deltam_cr_ct or pars::use_deltam_tm_ct)) CRASH("Need QED with mode 1 to include deltam");
      
      pars::can_change_pars=false;
      
      lprop::set_ins();
      qprop::set_ins();
      jqprop::set_ins();
      
      sigma::set_ins();
      pr_bil::set_ins();
      pr_meslep::set_ins();
    }
}
