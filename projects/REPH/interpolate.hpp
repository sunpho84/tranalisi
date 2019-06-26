#ifndef _INTERPOLATE_HPP
#define _INTERPOLATE_HPP

#include <permes_combo.hpp>

template <typename TV>
template <typename F,typename S>
auto permes_t<TV>::interpolate(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<S>& data)
{
  switch((m1.size()!=1)*2+(m2.size()!=1))
    {
    case 0:
      return interpolate0D(f,data);
      break;
    case 1:
      return interpolate1D(m2,f,data);
      break;
    case 2:
      return interpolate1D(m1,f,data);
      break;
    case 3:
    default:
      return interpolate2D(m1,m2,f,data);
      break;
    }
}

//! 2 masses interpolation
template <typename TV>
template <typename F,typename S>
auto permes_t<TV>::interpolate2D(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<S>& data)
{
  //cout<<"2D interpolation, masses "<<m1<<m2<<endl;
  
  //! Number of combinations
  const index_t indM({{"m1",m1.size()},{"m2",m2.size()}});
  
  const size_t nMcombos=indM.max();
  cout<<"nMCombos: "<<nMcombos<<endl;
  
  //! Type of the output
  typedef decltype(f(data[0])) R;
  
  //! Number of elements to be interpolated
  const size_t n=f(data[0]).size();
  
  //! Result
  R res(n,0.0);
  
  //! Independent variables
  vector<vector<double>> x(nMcombos,vector<double>(3));
  for(size_t iM=0;iM<nMcombos;iM++)
    {
      const vector<size_t> c=indM(iM);
      x[iM][0]=1;
      x[iM][1]=m1[c[0]];
      x[iM][2]=m2[c[1]];
    }
  
  for(size_t i=0;i<n;i++)
    {
      //! Type of an element of the output
      using Ri=typename remove_reference<decltype(res[0])>::type;
      
      //! Data for the fit
      plan_fit_data_t<Ri> plan_fit_data;
      
      for(size_t iM=0;iM<nMcombos;iM++)
	{
	  const auto t=make_tuple(x[iM],f(data[iM])[i]);
	  
	  plan_fit_data.push_back(t);
	}
      
      const auto coeffs=plan_fit(plan_fit_data);
      res[i]=plan_eval(coeffs,x[0]);
    }
  
  return res;
}

//! 1 mass interpolation
template <typename TV>
template <typename F,typename S>
auto permes_t<TV>::interpolate1D(const vector<double>& m,const F& f,const vector<S>& data)
{
  //cout<<"1D interpolation, masses "<<m<<endl;
  
  const size_t nM=m.size();
  
  //! Type of the output
  typedef decltype(f(data[0])) R;
  
  //! Number of elements to be interpolated
  const size_t n=f(data[0]).size();
  
  //! Result
  R res(n,0.0);
  
  //! Independent variables
  for(size_t i=0;i<n;i++)
    {
      R y(nM);
      //cout<<y<<endl;
      for(size_t iM=0;iM<nM;iM++)
	y[iM]=f(data[iM])[i];
      
      const auto coeffs=poly_fit(m,y,1);
      
      res[i]=poly_eval(coeffs,m[0]);
    }
  
  return res;
}

//! 0 mass interpolation
template <typename TV>
template <typename F,typename S>
auto permes_t<TV>::interpolate0D(const F& f,const vector<S>& data)
{
  //cout<<"0D interpolation"<<endl;
  return f(data[0]);
}

#endif
