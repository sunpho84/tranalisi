#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_COMMON
#include <common.hpp>

#include <set>

dboot_t read_boot(const raw_file_t &file)
{
  dboot_t out;
  for(size_t ib=0;ib<nboots;ib++) file.read(out[ib]);
  out.fill_ave_with_components_ave();
  return out;
}

// //! truncate to next most significative digit
double get_ml_max()
{
// //must be fixed because before commenting it was not getting data
// double ml_max=0;
// for(auto &data : ext_data)
//   ml_max=max(ml_max,dboot_t(data.aml/pars.fit_a[data.ib]/pars.fit_z[data.ib]).ave());
//   double l=floor(log(ml_max)/log(10));
//   double pad=pow(10,l);
//   ml_max=ceil(ml_max/pad)*pad;
//   return ml_max;
  return 0.05;
}

void init_common_IB(string ens_pars)
{
  set_njacks(15);
  
  raw_file_t file(ens_pars,"r");
  
  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  for(size_t ibeta=0;ibeta<nbeta;ibeta++)
    file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    {
      for(size_t iboot=0;iboot<nboots;iboot++)
	for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	  file.read(lat_par[input_an_id].ainv[ibeta][iboot]);
      for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].ainv[ibeta].fill_ave_with_components_ave();
    }
  file.expect({"r0","(GeV^-1)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].r0=read_boot(file);
  file.expect({"Zp","(1.90","1.95","2.10)"});
  for(int ihalf=0;ihalf<2;ihalf++)
    {
      //skip average of the half
      for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	file.read(dum);
      
      size_t batch=ninput_an/2;
      for(size_t input_an_id=batch*ihalf;input_an_id<batch*(ihalf+1);input_an_id++)
	{
	  for(size_t iboot=0;iboot<nboots;iboot++)
	    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
	      file.read(lat_par[input_an_id].Z[ibeta][iboot]);
	  for(size_t ibeta=0;ibeta<nbeta;ibeta++) lat_par[input_an_id].Z[ibeta].fill_ave_with_components_ave();
	}
    }
  file.expect({"Jackknife","numbers","(","0.0030(32),0.0040(32),","0.0050(32),","0.0040(24)","...)"});
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++)
    for(size_t iboot=0;iboot<nboots;iboot++)
      for(size_t iens=0;iens<nens_total;iens++)
	{
	  size_t ijack_plus_one;
	  file.read(ijack_plus_one);
	  jack_index[input_an_id][iens][iboot]=ijack_plus_one-1;
	}
  file.expect({"f0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].f0=read_boot(file);
  file.expect({"2*B0","(GeV)"});
  file.read(dum);
  for(size_t input_an_id=0;input_an_id<ninput_an;input_an_id++) lat_par[input_an_id].B0=read_boot(file)/2.0;
}

//! perform the analysis according to eq.28
ave_err_t eq_28_analysis(const dbvec_t &v)
{
  ave_err_t ae;
  double sigma=0;
  
  for(size_t i=0;i<v.size();i++)
    {
      double a=v[i].ave();
      double e=v[i].err();
      ae.ave()+=a;
      ae.err()+=sqr(a);
      sigma+=sqr(e);
    }
  ae.ave()/=v.size();
  ae.err()/=v.size();
  sigma/=v.size();
  ae.err()-=sqr(ae.ave());
  ae.err()=sqrt(fabs(ae.err())+sigma);
  
  return ae;
}

void cont_chir_fit_minimize
(const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,boot_fit_t &boot_fit,double apow,double zpow,
 const function<double(const vector<double> &p,const cont_chir_fit_pars_t &pars,double ml,double ms,double Maux,double ac,double L)> &cont_chir_ansatz,bool cov_flag)
{
  //set_printlevel(3);
  
  //set data
  for(size_t idata=0;idata<ext_data.size();idata++)
    boot_fit.add_point(//numerical data
		       [&ext_data,&pars,idata,apow,zpow]
		       (const vector<double> &p,int iel) //dimension 2
		       {return ext_data[idata].wfse[iel]*pow(pars.get_z(p,ext_data[idata].ib,iel),zpow)/pow(pars.get_a(p,ext_data[idata].ib,iel),apow);},
		       //ansatz
		       [idata,&pars,&ext_data,&cont_chir_ansatz]
		       (const vector<double> &p,int iel)
		       {
			 size_t ib=ext_data[idata].ib;
			 double ac=pars.get_a(p,ib,iel);
			 double zc=pars.get_z(p,ib,iel);
			 double ml=ext_data[idata].aml/ac/zc;
			 double ms=ext_data[idata].ams/ac/zc;
			 double Maux=ext_data[idata].aMaux[iel]/ac;
			 double L=ext_data[idata].L;
			 return cont_chir_ansatz(p,pars,ml,ms,Maux,ac,L);
		       },
		       //for covariance/error
		       dboot_t(ext_data[idata].wfse*pow(pars.ori_z[ext_data[idata].ib],zpow)/pow(pars.ori_a[ext_data[idata].ib],apow)),1/*correlate*/);
  
  //! fit
  boot_fit.fit(cov_flag);
  
  //print parameters
  pars.print_common_pars();
  pars.print_LEC_pars();
}

void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,const string &subtitle)
{
  //search max renormalized mass
  double ml_max=get_ml_max();
  
  //prepare plot
  grace_file_t fit_file(path);
  //fit_file.set_title("Continuum and chiral limit");
  fit_file.set_subtitle(subtitle);
  fit_file.set_xaxis_label("$$m_{light} (\\overline{MS},2 GeV) [GeV]");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  // grace::default_color_scheme={grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
  // grace::default_symbol_scheme={grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<pars.fit_a.size();ib++)
	{
	  //make the list of volumes
	  set<size_t> L_list;
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      L_list.insert(ext_data[idata].L);
	  
	  //loop over the list of volumes
	  for(auto &L : L_list)
	    {
	      fit_file.new_data_set();
	      //put data without fse to brown
	      if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	      if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
	      
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		  {
		    
		    fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		  }
	    }
	}
      
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
  fit_file.set_legend("physical point");
}

void plot_chir_fit_empty(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list)
{
  //search max renormalized mass
  double ml_max=0;
  for(auto &data : ext_data)
    ml_max=max(ml_max,dboot_t(data.aml/pars.fit_a[data.ib]/pars.fit_z[data.ib]).ave());
  ml_max*=1.1;
  
  //prepare plot
  grace_file_t fit_file(path);
  //fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("m\\s\\f{Times-Italic}l\\f{}\\N\\h{0.3}(GeV)");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_color_scheme={grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
  grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<pars.fit_a.size();ib++)
	{
	  //make the list of volumes
	  set<size_t> L_list;
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      L_list.insert(ext_data[idata].L);
	  
	  //loop over the list of volumes
	  for(auto &L : L_list)
	    {
	      fit_file.new_data_set();
	      //put data without fse to brown
	      if(without_with_fse==0) grace::default_symbol_fill_pattern=grace::EMPTY_SYMBOL;	      
	      if(without_with_fse==1)
		{
		  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
		  fit_file.set_legend(combine("$$\\beta=%s, L/a=%d",beta_list[ib].c_str(),L).c_str());
		}
	      
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		  {
		    
		    fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		  }
	    }
	}
      
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
  fit_file.set_legend("physical point");
}

void plot_chir_fit1(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list)
{
  //search max renormalized mass
  double ml_max=0;
  for(auto &data : ext_data)
    ml_max=max(ml_max,dboot_t(data.aml/pars.fit_a[data.ib]/pars.fit_z[data.ib]).ave());
  ml_max*=1.1;
  
  //prepare plot
  grace_file_t fit_file(path);
  //fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("m\\s\\f{Times-Italic}l\\f{}\\N\\h{0.3}(GeV)");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
  grace::default_color_scheme={grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
  grace::default_symbol_scheme={grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
  for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
    {
      for(size_t ib=0;ib<pars.fit_a.size();ib++)
	{
	  //make the list of volumes
	  set<size_t> L_list;
	  for(size_t idata=0;idata<ext_data.size();idata++)
	    if(ext_data[idata].ib==ib)
	      L_list.insert(ext_data[idata].L);
	  
	  //loop over the list of volumes
	  for(auto &L : L_list)
	    {
	      fit_file.new_data_set();
	      //put data without fse to brown
	      if(without_with_fse==0) fit_file.set_all_colors(grace::BROWN);
	      if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L/a=%d",beta_list[ib].c_str(),L).c_str());
	      
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		  {
		    
		    fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		  }
	    }
	}
      
      //put back colors for data with fse
      if(without_with_fse==0) fit_file.reset_cur_col();
    }
  //data of the continuum-chiral limit
  fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
  fit_file.set_legend("physical point");
}

void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,size_t univ_full_sub)
{
  //search max renormalized mass
  double ml_max=get_ml_max();
  
  //prepare plot
  grace_file_t fit_file(path);
  //fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$m_{light} (\\overline{MS},2 GeV) [GeV]");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  if(univ_full_sub==1)
    for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  if(univ_full_sub==1)
    fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  grace::default_color_scheme={grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
  grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
  if(univ_full_sub==0)
    {
      for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
	{
	  for(size_t ib=0;ib<pars.fit_a.size();ib++)
	    {
	      //make the list of volumes
	      set<size_t> L_list;
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib)
		  L_list.insert(ext_data[idata].L);
	  
	      //loop over the list of volumes
	      for(auto &L : L_list)	      
		{
		  fit_file.new_data_set();
		  //put data without fse to brown
		  if(without_with_fse==0)
		    {
		      grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
		      fit_file.set_all_colors(grace::BROWN);
		    }
		  if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
		        
		  for(size_t idata=0;idata<ext_data.size();idata++)
		    if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		      {
		    
			fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
			
		      }
		}
	    }
      
	  //put back colors for data with fse
	  if(without_with_fse==0)
	    {
	      grace::default_symbol_fill_pattern=grace::EMPTY_SYMBOL;
	      fit_file.reset_cur_col();
	    }
	}
    }
  else
    {
      for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
	{
	  fit_file.reset_cur_col();
	  for(size_t ib=0;ib<pars.fit_a.size();ib++)
	    {
	      //make the list of volumes
	      set<size_t> L_list;
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib)
		  L_list.insert(ext_data[idata].L);
	  
	      //loop over the list of volumes
	      for(auto &L : L_list)
		{
		  fit_file.new_data_set();
		  //put data without fse to brown
		  if(without_with_fse==0) grace::default_symbol_fill_pattern=grace::EMPTY_SYMBOL;
		  if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
	      
		  for(size_t idata=0;idata<ext_data.size();idata++)
		    if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		      {
		    
			fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		      }
		}
	    }
	  
	  //put back colors for data with fse
	  if(without_with_fse==0) grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
	      
	}
    }
  //data of the continuum-chiral limit
  if(univ_full_sub==1)
    {
      fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
      fit_file.set_legend("physical point");
    }
}

void plot_chir_fit(const string path,const vector<cont_chir_fit_data_t> &ext_data,const cont_chir_fit_pars_t &pars,
		   const function<double(double x,size_t ib)> &fun_line_per_beta,
		   const function<dboot_t(double x)> &fun_poly_cont_lin,
		   const function<dboot_t(size_t idata,bool without_with_fse,size_t ib)> &fun_data,
		   const dboot_t &ml_phys,const dboot_t &phys_res,const string &yaxis_label,const vector<string> &beta_list,size_t univ_full_sub,size_t FSE_flag)
{
  //search max renormalized mass
  double ml_max=get_ml_max();
  
  //prepare plot
  grace_file_t fit_file(path);
  //fit_file.set_title("Continuum and chiral limit");
  fit_file.set_xaxis_label("$$m_{light} (\\overline{MS},2 GeV) [GeV]");
  fit_file.set_yaxis_label(yaxis_label);
  fit_file.set_xaxis_max(ml_max);
  
  //band of the fit to individual beta
  if(univ_full_sub==1)
    for(size_t ib=0;ib<pars.fit_a.size();ib++) fit_file.write_line(bind(fun_line_per_beta,_1,ib),1e-6,ml_max);
  //band of the continuum limit
  if(univ_full_sub==1)
    fit_file.write_polygon(fun_poly_cont_lin,1e-6,ml_max);
  //data without and with fse
  if(FSE_flag==0)
    {
      grace::default_color_scheme={grace::RED,grace::BLUE,grace::GREEN4,grace::VIOLET};
      grace::default_symbol_scheme={grace::DIAMOND,grace::DIAMOND,grace::DIAMOND};
    }
  else
    {
      grace::default_color_scheme={grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET};
      grace::default_symbol_scheme={grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND};
    }
  if(univ_full_sub==0)
    {
      for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
	{
	  for(size_t ib=0;ib<pars.fit_a.size();ib++)
	    {
	      //make the list of volumes
	      set<size_t> L_list;
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib)
		  L_list.insert(ext_data[idata].L);
	      
	      //loop over the list of volumes
	      for(auto &L : L_list)
		{		
		  fit_file.new_data_set();
		  //put data without fse to brown
		  if(without_with_fse==0)
		    {
		      grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
		      fit_file.set_all_colors(grace::BROWN);
		    }
		  if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
		        
		  for(size_t idata=0;idata<ext_data.size();idata++)
		    if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		      {
		    
			fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		      }
		}
	    }
      
	  //put back colors for data with fse
	  if(without_with_fse==0)
	    {
	      grace::default_symbol_fill_pattern=grace::EMPTY_SYMBOL;
	      fit_file.reset_cur_col();
	    }
	}
    }
  else
    {
      for(int without_with_fse=0;without_with_fse<2;without_with_fse++)
	{
	  fit_file.reset_cur_col();
	  for(size_t ib=0;ib<pars.fit_a.size();ib++)
	    {
	      //make the list of volumes
	      set<size_t> L_list;
	      for(size_t idata=0;idata<ext_data.size();idata++)
		if(ext_data[idata].ib==ib)
		  L_list.insert(ext_data[idata].L);
	  
	      //loop over the list of volumes
	      for(auto &L : L_list)
		{
		  fit_file.new_data_set();
		  //put data without fse to brown
		  if(without_with_fse==0) grace::default_symbol_fill_pattern=grace::EMPTY_SYMBOL;
		  if(without_with_fse==1) fit_file.set_legend(combine("$$\\beta=%s, L=%d",beta_list[ib].c_str(),L).c_str());
	      
		  for(size_t idata=0;idata<ext_data.size();idata++)
		    if(ext_data[idata].ib==ib and ext_data[idata].L==L)
		      fit_file.write_ave_err(dboot_t(ext_data[idata].aml/pars.fit_z[ib]/pars.fit_a[ib]).ave(),fun_data(idata,without_with_fse,ib).ave_err());
		}
	    }
      
	  //put back colors for data with fse
	  if(without_with_fse==0) grace::default_symbol_fill_pattern=grace::FILLED_SYMBOL;
	      
	}
    }
  //data of the continuum-chiral limit
  if(univ_full_sub==1)
    {
      fit_file.write_ave_err(ml_phys.ave_err(),phys_res.ave_err());
      fit_file.set_legend("physical point");
    }
}
