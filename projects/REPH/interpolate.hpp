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
permes_t<dbvec_t> interpolate(const perens_t& ens,const string& mesName,const vector<double>& m1,const vector<double>& m2,const vector<permes_combo_t<djvec_t> >& data,const boot_init_t& j,const dboot_t& m1phys,const dboot_t& m2phys,const bool&forceInterp=false)
{
  //! Output
  permes_t<dbvec_t> out(ens,mesName);
  
  //! Path where to store or load
  const string dataPath=out.milledPath()+"/interp.dat";
  
  //! Load or interpolate
  const bool interpLoad=file_exists(dataPath) and not forceInterp;
  
  if(not interpLoad)
    {
      cout<<"Interpolating"<<endl;
      
#define INTERP(A)							\
      out.A=interpolate(m1,m2,[=](const auto& D){return D.A;},data,j,m1phys,m2phys)
      
      INTERP(E);
      INTERP(X);
      for(int i=0;i<2;i++)
	INTERP(ff[i]);
    }
#undef INTERP
  
  for(size_t iVA=0;iVA<2;iVA++)
    {
      out.quality[iVA].resize(ens.nDecKin,true);
      for(size_t iDecKin=0;iDecKin<ens.nDecKin;iDecKin++)
	for(auto& d : data)
	  out.quality[iVA][iDecKin]&=d.quality[iVA][iDecKin];
    }
  
  cout<<(interpLoad?"Reading":"Writing")<<" interpolated data"<<endl;
  //! Input or output file
  raw_file_t file(dataPath,interpLoad?"r":"w");
  
  for(auto& q: {&out.E,&out.X,&out.ff[0],&out.ff[1]})
    {
      if(interpLoad)
	{
	  q->resize(file.bin_read<size_t>());
	  file.bin_read(*q);
	}
      else
	{
	  file.bin_write(q->size());
	  file.bin_write(*q);
	}
    }
  
  return out;
}

//! Check if the range allows for an interpolation or consists of an extrapolation
void checkIntExtr(const vector<double>& amBare,const dboot_t& target)
{
  if(amBare.size()!=1)
    {
      //! Check if interpolation
      bool isInterpolating=true;
      
      const auto& amMinMax=minmax_element(amBare.begin(),amBare.end());
      const double& amMin=*amMinMax.first;
      const double& amMax=*amMinMax.second;
      const double range=amMax-amMin;
      
      size_t nExtrapolating=0;
      double maxExcess=0;
      
      for(size_t iboot=0;iboot<nboots;iboot++)
	{
	  const double thisBoot=target[iboot];
	  bool below=false,above=false;
	  for(const double& am : amBare)
	    {
	      if(am<thisBoot) below=true;
	      if(am>thisBoot) above=true;
	    }
	  const bool thisIsInterpolating=above and below;
	  
	  isInterpolating&=thisIsInterpolating;
	  if(not thisIsInterpolating)
	    {
	      nExtrapolating++;
	      
	      const double excess=
		(not below)?
		(amMin-thisBoot)
		:
		(thisBoot-amMax);
	      
	      maxExcess=max(maxExcess,excess);
	    }
	}
      
      if(nExtrapolating)
	cout<<"WARNING, extrapolating in: "<<nExtrapolating*100/nboots<<"% of the bootstraps, max extrapolation: "<<maxExcess*100/range<<"%"<<endl<<endl;
    }
}

//! Frontend to interpolate in mass
permes_t<dbvec_t>interpolate(const AllMesCombos& mesCombos,const meson_t& mesComposition,const perens_t& ens,const map<string,vector<double>>& am,const size_t& inputAn)
{
  //! Index of beta
  const size_t& iBeta=ens.iBeta;
  
  //! Name of the meson
  const string& mesName=get<0>(mesComposition)+"/interp/"+to_string(inputAn);
  
  //! S quark name
  const string& qS=get<1>(mesComposition);
  
  //! T quark name
  const string& qT=get<2>(mesComposition);
  
  //! Index of S quark in the physical list
  const size_t iMphysS=get<1>(quarkList[qS]);
  
  //! Index of T quark in the physical list
  const size_t iMphysT=get<1>(quarkList[qT]);
  
  //! Bare mass of S quark
  const dboot_t amSbare=lat_par[inputAn].amBare(iMphysS,iBeta);
  
  //! Bare mass of T quark
  const dboot_t amTbare=lat_par[inputAn].amBare(iMphysT,iBeta);
  
  for(auto& i : {make_tuple("S",qS,amSbare),make_tuple("T",qT,amTbare)})
    {
      const auto& tag=get<0>(i);
      const vector<double> amList=am.at(get<1>(i));
      const dboot_t& amBare=get<2>(i);
      
      cout<<"am"<<tag<<"bare: "<<smart_print(amBare);
      for(auto& ami : amList)
	cout<<" "<<ami;
      cout<<endl;
      checkIntExtr(amList,amBare);
    }
  
  return interpolate(ens,mesName,am.at(qS),am.at(qT),mesCombos,jack_index[inputAn][ens.iUlt],amSbare,amTbare);
}

#endif
