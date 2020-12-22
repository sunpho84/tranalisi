#ifndef _NODE_HPP
#define _NODE_HPP

#include <vector>
#include <ostream>
#include <set>
#include <sstream>

#include <tranalisi.hpp>

using namespace std;

/// Generic node
struct node_t
{
  /// Virtual evaluator
  virtual double eval(const vector<double>& pars,const double& x) const = 0;
  
  /// Virtual destructor
  virtual ~node_t() = 0;
  
  /// Print to an ostream
  virtual ostream& print(ostream&) const = 0;
  
  /// Get parameters
  virtual set<size_t> getParIds() const = 0;
  
  /// Evaluate over the jackknife
  djack_t jackEval(const djvec_t& pars,const double& x) const
  {
    /// Result
    djack_t y;
    
    /// Number of parameters
    const size_t nPars=pars.size();
    
    /// Parameters
    vector<double> parsEl(nPars);
    
    for(size_t ijack=0;ijack<=njacks;ijack++)
      {
	for(size_t iPar=0;iPar<nPars;iPar++)
	  parsEl[iPar]=pars[iPar][ijack];
	
	y[ijack]=this->eval(parsEl,x);
      }
    
    return y;
  }
};

/// Virtual destructor
inline node_t::~node_t()
{
}

/////////////////////////////////////////////////////////////////

/// Variable
struct variable_node_t : public node_t
{
  /// Creates
  variable_node_t()
  {
  }
  
  /// Destruct
  ~variable_node_t()
  {
  }
  
  /// Evaluate
  double eval(const vector<double>& pars,const double& x) const
  {
    return x;
  }
  
  /// Print to an ostream
  ostream& print(ostream& os) const
  {
    return os<<"x";
  }
  
  /// Returns the parameters
  set<size_t> getParIds() const
  {
    return {};
  };
};

/////////////////////////////////////////////////////////////////

/// Parameter
struct parameter_node_t : public node_t
{
  /// Reference in the par list
  const size_t id;
  
  /// Creates taking the id
  parameter_node_t(size_t id) : id(id)
  {
  }
  
  /// Destruct
  ~parameter_node_t()
  {
  }
  
  /// Evaluate
  double eval(const vector<double>& pars,const double& x) const
  {
    return pars[id];
  }
  
  /// Print to an ostream
  ostream& print(ostream& os) const
  {
    return os<<"["<<id<<"]";
  }
  
  /// Returns the parameters
  set<size_t> getParIds() const
  {
    return {id};
  };
};

/////////////////////////////////////////////////////////////////

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
  double eval(const vector<double>& pars,const double& x) const
  {
    return val;
  }
  
  /// Print to an ostream
  ostream& print(ostream& os) const
  {
    return os<<val;
  }
  
  /// Returns the parameters
  set<size_t> getParIds() const
  {
    return {};
  };
};

/////////////////////////////////////////////////////////////////

/// Unary expression
struct uexp_node_t : public node_t
{
  /// Function to be called
  double (*fun)(const double);
  
  /// Full name
  const string name;
  
  /// Argument of the function
  const node_t* arg;
  
  /// Creates taking the function and the argument
  uexp_node_t(double (*fun)(const double),const string& name,const node_t* arg) : fun(fun),name(name),arg(arg)
  {
  }
  
  /// Destruct subnodes
  ~uexp_node_t()
  {
    delete arg;
  }
  
  /// Evaluate
  double eval(const vector<double>& pars,const double& x) const
  {
    return fun(arg->eval(pars,x));
  }
  
  /// Print to an ostream
  ostream& print(ostream& os) const
  {
    ostringstream oss;
    arg->print(oss);
    
    return os<<combine(name.c_str(),oss.str().c_str());
  }
  
  /// Returns the parameters
  set<size_t> getParIds() const
  {
    return arg->getParIds();
  };
};

/////////////////////////////////////////////////////////////////

/// Binary expression
struct bexp_node_t : public node_t
{
  /// Function to be called
  double (*fun)(const double,const double);
  
  /// Full name of the function
  const string name;
  
  /// First argument of the function
  const node_t* arg1;
  
  /// Second argument of the function
  const node_t* arg2;
  
  /// Creates taking the function and the argument
  bexp_node_t(double (*fun)(const double,const double),const string& name,const node_t* arg1,const node_t* arg2) : fun(fun),name(name),arg1(arg1),arg2(arg2)
  {
  }
  
  /// Destruct subnodes
  ~bexp_node_t()
  {
    delete arg1;
    delete arg2;
  }
  
  /// Evaluate
  double eval(const vector<double>& pars,const double& x) const
  {
    return fun(arg1->eval(pars,x),arg2->eval(pars,x));
  }
  
  /// Print to an ostream
  ostream& print(ostream& os) const
  {
    ostringstream oss1,oss2;
    arg1->print(oss1);
    arg2->print(oss2);
    
    return os<<combine(name.c_str(),oss1.str().c_str(),oss2.str().c_str());
  }
  
  /// Returns the parameters
  set<size_t> getParIds() const
  {
    set<size_t> res=arg1->getParIds();
    for(auto& i : arg2->getParIds())
      res.insert(i);
    
    return res;
  };
};

/// Print to an ostream
inline ostream& operator<<(ostream& os,const node_t* node)
{
  return node->print(os);
}

#endif
