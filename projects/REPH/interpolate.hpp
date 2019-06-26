#ifndef _INTERPOLATE_HPP
#define _INTERPOLATE_HPP

#include <permes_combo.hpp>

//! 0 1 or 2D version frontend
template <typename F>
dbvec_t interpolate(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<permes_combo_t<djvec_t>>& data,const boot_init_t& j,const dboot_t& m1phys,const dboot_t& m2phys)
{
  switch((m1.size()!=1)*2+(m2.size()!=1))
    {
    case 0:
      return interpolate0D(f,data,j);
      break;
    case 1:
      return interpolate1D(m2,f,data,j,m2phys);
      break;
    case 2:
      return interpolate1D(m1,f,data,j,m1phys);
      break;
    case 3:
    default:
      return interpolate2D(m1,m2,f,data,j,m1phys,m2phys);
      break;
    }
}

//! 2 masses interpolation
template <typename F>
dbvec_t interpolate2D(const vector<double>& m1,const vector<double>& m2,const F& f,const vector<permes_combo_t<djvec_t>>& data,const boot_init_t& j,const dboot_t& m1phys,const dboot_t& m2phys)
{
  //! Number of combinations
  const index_t indM({{"m1",m1.size()},{"m2",m2.size()}});
  
  const size_t nMcombos=indM.max();
  
  //! Number of elements to be interpolated
  const size_t n=f(data[0]).size();
  
  //! Result
  dbvec_t res(n);
  
  //! Independent variables
  vector<vector<double>> x(nMcombos,vector<double>(3));
  for(size_t iM=0;iM<nMcombos;iM++)
    {
      const vector<size_t> c=indM(iM);
      x[iM][0]=1;
      x[iM][1]=m1[c[0]];
      x[iM][2]=m2[c[1]];
    }
  
  dbvec_t xphys(3);
  xphys[0]=1;
  xphys[1]=m1phys;
  xphys[2]=m2phys;
  
  for(size_t i=0;i<n;i++)
    {
      //! Data for the fit
      plan_fit_data_t<dboot_t> plan_fit_data;
      
      for(size_t iM=0;iM<nMcombos;iM++)
	{
	  const auto t=make_tuple(x[iM],dboot_t(j,f(data[iM])[i]));
	  
	  plan_fit_data.push_back(t);
	}
      
      const auto coeffs=plan_fit(plan_fit_data);
      res[i]=plan_eval(coeffs,xphys);
    }
  
  return res;
}

//! 1 mass interpolation
template <typename F,typename S>
dbvec_t interpolate1D(const vector<double>& m,const F& f,const vector<S>& data,const boot_init_t& j,const dboot_t& mphys)
{
  const size_t nM=m.size();
  
  //! Number of elements to be interpolated
  const size_t n=f(data[0]).size();
  
  //! Result
  dbvec_t res(n);
  
  //! Independent variables
  for(size_t i=0;i<n;i++)
    {
      dbvec_t y(nM);
      
      for(size_t iM=0;iM<nM;iM++)
	y[iM].fill_from_jack(j,f(data[iM])[i]);
      
      const auto coeffs=poly_fit(m,y,1);
      
      res[i]=poly_eval(coeffs,mphys);
    }
  
  return res;
}

//! 0 mass interpolation
template <typename F,typename S>
dbvec_t interpolate0D(const F& f,const vector<S>& data,const boot_init_t& j)
{
  return dbvec_t(j,f(data[0]));
}

//! Interpolate from a list
permes_t<dbvec_t> interpolate(const perens_t& ens,const string& mesName,const vector<double>& m1,const vector<double>& m2,const vector<permes_combo_t<djvec_t> >& data,const boot_init_t& j,const dboot_t& m1phys,const dboot_t& m2phys)
{
  permes_t<dbvec_t> out(ens,mesName);
  
#define INTERP(A)						\
  out.A=interpolate(m1,m2,[](const auto& D){return D.A;},data,j,m1phys,m2phys)
  
  INTERP(E);
  INTERP(X);
  INTERP(ff[0]);
  INTERP(ff[1]);
  
#undef INTERP
  
  return out;
}

#endif
