#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_ANALYSIS
 #include <MOM2/analysis.hpp>

#include <MOM2/Zq.hpp>

#define ASSERT_COMPATIBLE_MEMBER(MEMBER)			\
  if(data(in1,ASSERT_PRESENT).MEMBER!=				\
     data(in2,ASSERT_PRESENT).MEMBER)				\
    CRASH("Impossible to combine, different member " #MEMBER)

void assert_compatible(const string in1,const string in2)
{
  for(size_t mu=0;mu<NDIM;mu++) ASSERT_COMPATIBLE_MEMBER(L[mu]);
  ASSERT_COMPATIBLE_MEMBER(nm);
  ASSERT_COMPATIBLE_MEMBER(nr);
  ASSERT_COMPATIBLE_MEMBER(linmoms);
  ASSERT_COMPATIBLE_MEMBER(bilmoms);
}

void make_Z_QED_absolute()
{
  CRASH("Are you sure that Z in QCD do not contain an offset?");
  
  needs_to_read_Z();
  if(not pars::Z_QED_are_relative) CRASH("Need to have relative QED before");
  
  for(auto &path : pars::ens)
    data(path,ASSERT_PRESENT)
      .make_Z_QED_absolute();
  
  pars::Z_QED_are_relative=false;
  invalidate_ingredients();
}

void plot_all_Z(const string &suffix)
{
  needs_to_read_Z();
  
  for(auto &path : pars::ens)
    data(path,ASSERT_PRESENT)
      .plot_Z(suffix);
}

void print_all_Z(const string &path)
{
  needs_to_read_Z();
  
  ofstream file(path);
  if(not file.good()) CRASH("Opening %s",path.c_str());
  
  for(auto &path : pars::ens)
    {
      file<<path<<endl;
      
      data(path,ASSERT_PRESENT)
	.print_Z(file);
    }
}

void average_Z(const string out,const string in1,const string in2)
{
  combine_Z(out,in1,in2,"average",[](djvec_t& out,const djvec_t& in1,const djvec_t& in2){out=(in1+in2)/2.0;});
}

void average_ingredients(const string out_name,const string in1,const string in2)
{
  invalidate_Z();
  
  needs_to_read_ingredients();
  
  assert_compatible(in1,in2);
  
  cout<<"Averaging ingredients for "<<in1<<" and "<<in2<<" into "<<out_name<<endl;
  
  perens_t temp=data(in1,ASSERT_PRESENT);
  
  for(auto &p : temp.get_all_ingredients({&data(in2,ASSERT_PRESENT)}))
    {
      cout<<" "<<p.tag<<endl;
      
      djvec_t &out=*p.out;
      const djvec_t &in=*p.in.front();
      
      out+=in;
      out/=2.0;
    }
  auto &dout=temp;
  const auto &din1=data(in1,PRESENCE_NOT_NEEDED);
  const auto &din2=data(in2,PRESENCE_NOT_NEEDED);
  dout.meson_mass_sea=(din1.meson_mass_sea+din2.meson_mass_sea)/2.0;
  for(size_t im=0;im<dout.nm;im++)
    for(size_t r=0;r<dout.nr;r++)
      {
	size_t imr=dout.im_r_ind({im,r});
	cout<<"deltam_cr[m="<<im<<",r="<<r<<"]: "<<
	  smart_print(dout.deltam_cr[imr].ave_err())<<" = [ "<<
	  smart_print(din1.deltam_cr[imr].ave_err())<<" + "<<
	  smart_print(din2.deltam_cr[imr].ave_err())<<" )/2.0"<<endl;
	cout<<"deltam_tm[m="<<im<<",r="<<r<<"]: "<<
	  smart_print(dout.deltam_tm[imr].ave_err())<<" = [ "<<
	  smart_print(din1.deltam_tm[imr].ave_err())<<" + "<<
	  smart_print(din2.deltam_tm[imr].ave_err())<<" ]/2.0"<<endl;
      }
  
  for(size_t im1=0;im1<dout.nm;im1++)
    for(size_t im2=0;im2<dout.nm;im2++)
      {
	size_t i=dout.im_im_ind({im1,im2});
	cout<<"meson_mass[m1="<<im1<<",m2="<<im2<<"]: "<<
	  smart_print(dout.meson_mass[i].ave_err())<<" = ["<<
	  smart_print(din1.meson_mass[i].ave_err())<<" + "<<
	  smart_print(din2.meson_mass[i].ave_err())<<"] /2.0"<<endl;
      }
  cout<<"meson_mass_sea: "<<
    dout.meson_mass_sea.ave_err()<<" = "<<din1.meson_mass_sea.ave_err()<<" + "<<din2.meson_mass_sea.ave_err()<<endl;
  
  //remove from the list
  for(auto in : {in1,in2})
    data_erase(in);
  
  //add to the list
  pars::ens.push_back(out_name);
  data(out_name,PRESENCE_NOT_NEEDED)=temp;
  data(out_name,PRESENCE_NOT_NEEDED).dir_path=out_name;
}

void ratio_Z_minus_one(const string out,const string in1,const string in2)
{
  combine_Z(out,in1,in2,"ratio minus one",[](djvec_t& out,const djvec_t& in1,const djvec_t& in2){out=in1/in2-1.0;});
}

void subtract_Z(const string out,const string in1,const string in2)
{
  combine_Z(out,in1,in2,"ratio minus one",[](djvec_t& out,const djvec_t& in1,const djvec_t& in2){out=in1-in2;});
}

void sea_chir_extrap(const string out_name,const vector<string> &ens_list)
{
  needs_to_read_Z();
  
  invalidate_ingredients();
  
  if(ens_list.size()<2) CRASH("Need at least 2 ensembles, %zu passed",ens_list.size());
  cout<<"Chirally extrapolating "<<out_name<<" in the sea mass, list of ensembles:"<<endl;
  
  //take x and ensembles
  vector<double> x;
  vector<const perens_t*> in;
  for(auto &en_name : ens_list)
    {
      //check and take the data, add it to the list
      assert_compatible(ens_list.front(),en_name);
      const perens_t& en=data(en_name,ASSERT_PRESENT);
      in.push_back(&en);
      
      //compute the x
      double t;
      if(pars::chir_extr_method==chir_extr::MQUARK) t=en.am[en.im_sea];
      else                                          t=sqr(en.meson_mass_sea.ave());
      x.push_back(t);
      
      //print the name
      cout<<" "<<en_name<<endl;
    }
  
  auto xminmax=minmax_element(x.begin(),x.end());
  double xmin=*xminmax.first;
  xmin=0.0;
  double xmax=*xminmax.second*1.1;
  
  //prepare output
  perens_t out=*in.front();
  out.meson_mass_sea=0.0;
  out.dir_path=out_name;
  
  for(auto &v : out.get_all_Ztasks(in))
    {
      cout<<" "<<v.tag<<endl;
      
      for(size_t icombo=0;icombo<v.out->size();icombo++)
	{
	  djvec_t y(ens_list.size());
	  for(size_t iens=0;iens<v.in.size();iens++) y[iens]=(*(v.in[iens]))[icombo];
	  
	  string plot_path="";
	  if(icombo==100) plot_path=out_name+"/plots/sea_chirextr_"+v.tag+"_"+v.ind.descr(icombo)+".xmg";
	  djvec_t coeffs=poly_fit(x,y,1);
	  if(std::isnan(coeffs[0][0])) coeffs=0.0;
	  (*v.out)[icombo]=coeffs[0];
	  
	  if(plot_path!="")
	    {
	      grace_file_t plot(plot_path);
	      write_fit_plot(plot,xmin,xmax,bind(poly_eval<djvec_t>,coeffs,_1),x,y);
	      plot.set_title(v.tag+", "+v.ind.descr(icombo));
	      plot.write_ave_err(0,coeffs[0].ave_err());
	    }
	}
    }
  
  //remove from the list
  for(auto in : ens_list)
    data_erase(in);
  
  //add to the list
  pars::ens.push_back(out_name);
  data(out_name,PRESENCE_NOT_NEEDED)=out;
  data(out_name,PRESENCE_NOT_NEEDED).dir_path=out_name;
}

auto assert_ens_present(const string &key)
{
  auto f=_data.find(key);
  
  if(f==_data.end())
    {
      cout<<"Available keys: "<<endl;
      for(auto d : _data) cout<<" "<<d.first<<endl;
      CRASH("Unable to find the key %s",key.c_str());
    }
  
  return f;
}

void add_ens(const string &name)
{
  freeze_pars();
  
  pars::ens.push_back(name);
  
  perens_t &ens=data(name,PRESENCE_NOT_NEEDED);
  
  ens.read_pars(name)
    .set_pars_for_scratch()
    .set_indices();
  
  //allocate if asked to do immediately
  if(pars::allocate_immediately) ens.allocate();
}

void data_erase(const string &key)
{
  auto f=_data.find(key);
  auto e=find(pars::ens.begin(),pars::ens.end(),key);
  
  if(f!=_data.end() and e!=pars::ens.end())
    {
      _data.erase(key);
      pars::ens.erase(e);
    }
  else
    cout<<"Warning, key \""<<key<<"\" not present"<<endl;
}

perens_t& data(const string &key,const bool assert_present_flag)
{
  if(assert_present_flag)
    return assert_ens_present(key)->second;
  else
    return _data[key];
}

void list_ensembles()
{
  cout<<"Ensembles:"<<endl;
   for(auto &path : pars::ens)
     cout<<" "<<path<<endl;
}

void print_discr()
{
  for(auto &path : pars::ens) data(path,ASSERT_PRESENT).print_discr();
}

void compute_or_load_all_ingredients()
{
  cout<<"Going to rotate propagators: "<<(pars::twisted_run and pars::phys_basis)<<endl;
  
  for(auto &name : pars::ens)
    {
      perens_t &ens=data(name,ASSERT_PRESENT);
      
      //allocate if not asked to do immediately at definition
      if(not pars::allocate_immediately) ens.allocate();
      
      if(pars::report_mPCAC) ens.get_mPCAC();
      
      ens.read_or_compute_ingredients();
      if(pars::use_deltam_cr_ct or pars::use_deltam_tm_ct) ens.get_deltam();
      ens.get_meson_mass();
      
      if(pars::average_equiv_momenta_immediately) ens.average_equiv_momenta();
      if(pars::average_r_immediately) ens.average_r();
    }
  
  validate_ingredients();
}

void combined_sea_chir_extrap(const vector<comb_extr_t> &list)
{
  needs_to_read_Z();
  
  invalidate_ingredients();
  
  //check that nmom==1
  for(auto &l : list)
    for(auto &name : get<2>(l))
      {
	const size_t nmoms=data(name,ASSERT_PRESENT).all_moms.size();
	if(nmoms!=1) CRASH("%s has more than 1 momenta, %zu",name.c_str(),nmoms);
      }
  
  //get and check ngroups
  const size_t ngroups=list.size();
  if(ngroups<2) CRASH("Need at least 2 set of ensembles, %zu passed",ngroups);
  
  //input data per each ensemble
  struct input_per_ens_t
  {
    double m2;    //the ensemble par
    djack_t in;   //the input
  };
  
  //! data per group
  struct per_group_t
  {
    djack_t* out;                     //reference to the output
    double a;                         //the group par
    vector<input_per_ens_t> ens_list; //the input
  };
  
  //! data to be fitted, per quantity
  struct per_quantity_t
  {
    vector<per_group_t> in_out;       //one per group, the input and output
    per_quantity_t(const size_t ngroups) : in_out(ngroups) {}
  };
  
  //data to be fitted, all quantities sorted by tag and pars
  map<string,per_quantity_t> fit_data;
  for(size_t igroup=0;igroup<ngroups;igroup++)
    {
      auto& l=list[igroup];
      const string &name_out=get<0>(l);
      const double &a=get<1>(l);
      const vector<string> &names_in=get<2>(l);
      
      //prepare output
      cout<<"Out: "<<name_out<<endl;
      pars::ens.push_back(name_out);
      data(name_out,PRESENCE_NOT_NEEDED)=data(names_in.front(),ASSERT_PRESENT);
      perens_t& data_out=data(name_out,ASSERT_PRESENT);
      data_out.dir_path=name_out;
      
      //get data list, check compatiblity and copy first in list
      vector<const perens_t*> data_in;
      for(size_t i=0;i<names_in.size();i++)
	{
	  assert_compatible(name_out,names_in[i]);
	  data_in.push_back(&data(names_in[i],ASSERT_PRESENT));
	  
	  cout<<" "<<names_in[i]<<endl;
	}
      
      //transpose data
      for(auto &p : data_out.get_all_Ztasks(data_in))
	{
	  //prepare output
	  for(size_t i=0;i<p.ind.max();i++)
	    {
	      //prepare the list referring to a "static"
	      const string tag=p.tag+"_"+p.ind.descr(i);
	      auto _d=fit_data.find(tag);
	      if(_d==fit_data.end()) fit_data.emplace(tag,per_quantity_t(ngroups));
	      auto& d=fit_data.find(tag)->second;
	      per_group_t& in_out=d.in_out[igroup];
	      in_out.out=&((*p.out)[i]);
	      
	      //prepare input
	      const size_t nens=p.in.size();
	      in_out.ens_list.resize(nens);
	      in_out.a=a;
	      for(size_t iens=0;iens<nens;iens++)
		{
		  //get mass
		  const djvec_t &y=(*(p.in[iens]));
		  const perens_t& en=data(names_in[iens],ASSERT_PRESENT);
		  double m2=sqr(en.meson_mass_sea.ave())/sqr(a);
		  
		  //prepare data per ens
		  input_per_ens_t& input=in_out.ens_list[iens];
		  input.m2=m2;
		  input.in=y[i];
		}
	    }
	}
    }
  
  //extrapolate all quantities
  for(auto &m : fit_data)
    {
      const string &tag=m.first;
      const per_quantity_t &per_quantity=m.second;
      
      cout<<"Extrapolating "<<tag<<" to chiral limit in the sea"<<endl;
      
      //fit and pars
      jack_fit_t jack_fit;
      djack_t mslope=0.0;
      djack_t mslope_a2dep=0.0;
      
      //output plot
      grace_file_t plot("comb_sea_chir/plots/comb_sea_chir_extrap_"+tag+".xmg");
      {
	using namespace grace;
	auto color_scheme={RED,RED,BLUE,BLUE,GREEN4,GREEN4,VIOLET,VIOLET};
	plot.set_color_scheme(color_scheme);
	plot.set_line_color_scheme(color_scheme);
	plot.set_symbol_scheme({SQUARE,DIAMOND,SQUARE,DIAMOND,SQUARE,DIAMOND,SQUARE,DIAMOND});
      }
      
      double m2_max=0.0;
      
      //define parameters
      vector<size_t> iextr(ngroups);
      size_t imslope=jack_fit.add_fit_par(mslope,"mslope",0.0,0.1);
      //size_t imslope_a2dep=jack_fit.add_fit_par(mslope_a2dep,"mslope_a2dep",0.0,0.1);
      //jack_fit.fix_par(imslope_a2dep);
      for(size_t igroup=0;igroup<ngroups;igroup++)
	{
	  const per_group_t per_group=per_quantity.in_out[igroup];
	  const vector<input_per_ens_t> &ens_list=per_group.ens_list;
	  //const double a=per_group.a;
	  
	  iextr[igroup]=jack_fit.add_fit_par(*per_group.out,"extr_"+to_string(igroup),0.0,0.1);
	  
	  plot.new_data_set();
	  
	  for(size_t iens=0;iens<ens_list.size();iens++)
	    {
	      const djack_t& fit_in=ens_list[iens].in;
	      const double& m2=ens_list[iens].m2;
	      m2_max=max(m2,m2_max);
	      jack_fit.add_point(//numerical data
				 [fit_in]
				 (const vector<double> &p,int iel)
				 {
				   return fit_in[iel];
				 },
				 //ansatz
				 [=]
				 (const vector<double> &p,int iel)
				 {
				   return p[iextr[igroup]]+m2*(p[imslope]// +sqr(a)*p[imslope_a2dep]
							       );
				 },
				 //for covariance/error
				 fit_in.err());
	      
	      plot.write_ave_err(m2,fit_in.ave_err());
	    }
	}
      
      //fit
      jack_fit.fit();
      
      //line of the fit
      for(size_t igroup=0;igroup<ngroups;igroup++)
	{
	  const djack_t& c0=*per_quantity.in_out[igroup].out;
	  const double a=per_quantity.in_out[igroup].a;
	  plot.write_line([&](const double x)
			  {
			    const djack_t out=c0+x*(mslope+sqr(a)*mslope_a2dep);
			    return out.ave();
			  },
			  0.0,m2_max);
	  plot.set_legend(get<0>(list[igroup]));
	}
      
      //extrapolated results
      plot.new_data_set();
      for(size_t igroup=0;igroup<ngroups;igroup++)
	{
	  const per_group_t per_group=per_quantity.in_out[igroup];
	  const djack_t &out=*per_group.out;
	  
	  plot.write_ave_err(0.0,out.ave_err());
	}
    }
  
  //remove all extrapolated ensembles
  for(auto &l : list)
    for(auto &name : get<2>(l))
      data_erase(name);
}
