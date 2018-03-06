#ifndef _DRIVER_HPP
#define _DRIVER_HPP

#include <fstream>
#include <iostream>

#include <ingredients.hpp>
#include <macros.hpp>
#include <tools.hpp>

using namespace std;

class driver_t
{
public:
  
  void *scanner;
  
  ifstream fin;
  
  driver_t(const string &path);
  
  void parse();
  
  //destructor
  ~driver_t()
  {
    destroy_scanner();
  }
  
  //! number of momenta between each print
  size_t print_each_mom;
  //! gauge action
  gaz::type_t act;
  //! beta
  double beta;
  
  map<string,ingredients_t> ing;
  string suff;
  
protected:
  void init_scanner();
  void destroy_scanner();
};

int parser_parse(driver_t *driver);

#endif

