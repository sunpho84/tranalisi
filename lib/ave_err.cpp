#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <ave_err.hpp>
#include <grace.hpp>
#include <oper.hpp>

void vec_ave_err_t::write(const string &path)
{grace_file_t(path)<<(*this);}

