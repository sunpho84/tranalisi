#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <ave_err.hpp>
#include <grace.hpp>
#include <oper.hpp>

void vec_ave_err_t::write(const string &path) const
{
  grace_file_t gr(path);
  gr.set_settype(grace::XYDY);
  for(size_t it=0;it<this->size();it++)
    if(!std::isnan((*this)[it].err()))
      gr<<it<<" "<<(*this)[it]<<endl;
  gr.new_data_set();
}

