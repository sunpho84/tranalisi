#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <grace.hpp>

void vec_ave_err_t::write(const string &path)
{
  grace_file_t file(path);
  file<<(*this);
  file.new_set();
}

