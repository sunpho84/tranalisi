#ifndef _NODE_HPP
#define _NODE_HPP

#include <vector>

using namespace std;

extern vector<double> pars;

/// Generic node
struct node_t
{
  /// Virtual evaluator
  virtual double eval() const = 0;
  
  /// Virtual destructor
  virtual ~node_t() = 0;
};

/// Virtual destructor
inline node_t::~node_t()
{
}

/// Parameter
struct parameter_node_t : public node_t
{
  /// Reference in the par list
  const int id;
  
  /// Creates taking the id
  parameter_node_t(int id) : id(id)
  {
  }
  
  /// Destruct
  ~parameter_node_t()
  {
  }
  
  /// Evaluate
  double eval() const
  {
    return pars[id];
  }
};

/// Real number
struct real_node_t : public node_t
{
  /// Value
  const double val;
  
  /// Creates taking the value
  real_node_t(double val) : val(val)
  {
  }
  
  /// Destruct
  ~real_node_t()
  {
  }
  
  /// Evaluate
  double eval() const
  {
    return val;
  }
};

/// Unary expression
struct uexp_node_t : public node_t
{
  /// Function to be called
  double (*fun)(const double);
  
  /// Argument of the function
  const node_t* arg;
  
  /// Creates taking the function and the argument
  uexp_node_t(double (*fun)(const double),const node_t* arg) : fun(fun),arg(arg)
  {
  }
  
  /// Destruct subnodes
  ~uexp_node_t()
  {
    delete arg;
  }
  
  /// Evaluate
  double eval() const
  {
    return fun(arg->eval());
  }
};

/// Binary expression
struct bexp_node_t : public node_t
{
  /// Function to be called
  double (*fun)(const double,const double);
  
  /// First argument of the function
  const node_t* arg1;
  
  /// Second argument of the function
  const node_t* arg2;
  
  /// Creates taking the function and the argument
  bexp_node_t(double (*fun)(const double,const double),const node_t* arg1,const node_t* arg2) : fun(fun),arg1(arg1),arg2(arg2)
  {
  }
  
  /// Destruct subnodes
  ~bexp_node_t()
  {
    delete arg1;
    delete arg2;
  }
  
  /// Evaluate
  double eval() const
  {
    return fun(arg1->eval(),arg2->eval());
  }
};

#endif
