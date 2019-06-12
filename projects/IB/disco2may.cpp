#include <tranalisi.hpp>

#include <set>

int T,L;
int TH;

//parametrization of the diagras
// enum {EU1,EU2,EU4,EU5,EU6};
const int ndiag=5;
const int diag[ndiag]={1,2,4,5,6};

vector<dcompl_t> read_vector(const string &path,int n,int nskip=0)
{
  raw_file_t data(path,"r");
  
  //cout<<"Reading from "<<path<<" "<<n<<" lines after skipping "<<nskip<<" words"<<endl;
  for(int iskip=0;iskip<nskip;iskip++)
    data.skip_line();
  
  vector<dcompl_t> d;
  
  for(int i=0;i<n;i++)
    {
      double re=data.read<double>();
      double im=data.read<double>();
      
      d.emplace_back(re,im);
    }
  
  return d;
}

template <typename T>
T average(const vector<T>& v)
{
  return accumulate(v.begin(),v.end(),dcompl_t{0,0})/(double)v.size();
}

///Source number
int Ndd(int div_nhits, int idiag){
  int b=(idiag<3)?div_nhits:div_nhits*(div_nhits-1);
  return b;
}
///compute divisors
vector<int> divisors(const int n)
{
  vector<int> b;
  
  for(int idiv=1;idiv<n+1;idiv++)
    if(n%idiv==0)
      b.push_back(idiv);
  
  return b;
}

vector<size_t> determine_list_of_confs(const range_t& conf_range)
{
  //determine the list of confs
  const vector<size_t> confs_p=get_existing_paths_in_range("out/%04d/mes_contr_Pion",conf_range,SILENT);
  const vector<size_t> confs_1=get_existing_paths_in_range("out/%04d/EU1_stoch",conf_range,SILENT);
  vector<size_t> confs;
  
  for(auto &_p : confs_p)
    for(auto & _1 : confs_1)
      if(_p==_1)
	  confs.push_back(_p);
  
  return confs;
}

template <typename F>
void slice_plot(const vector<int>& dim1,const int dim1_min,const vector<int>& dim2,const int dim2_min,const size_t idiag,const vector<vector<double>>& slopes,const vector<vector<double>>& errors,const index_t& ind,int idim1,int idim2,const F& fun,const double pow1,const double pow2)
{
  const char* name=ind.name(idim2).c_str();
  grace_file_t plot_ave(combine("plots/fit_results/EU%d_average_vs_%s",diag[idiag],name));
  grace_file_t plot_err(combine("plots/fit_results/EU%d_errors_vs_%s",diag[idiag],name));
  plot_ave.set_title(combine("Slopes vs %s",name));
  plot_err.set_title(combine("Errors vs %s",name));
  bool ytitle[5]={1,1,0,0,0};
  
  if(ytitle[idiag]==1)
    {
      plot_ave.set_yaxis_label("slopes");
      plot_err.set_yaxis_label("errors");
    }
  else
    {
      plot_ave.set_yaxis_label("slopes (a*M)");
      plot_err.set_yaxis_label("errors (a*M)");
    }
  for(grace_file_t* _g : {&plot_ave,&plot_err})
    {
      grace_file_t& g=(*_g);
      
      using namespace grace;
      
      g.set_subtitle(combine("EU%d",diag[idiag]));
      g.set_xaxis_logscale();
      g.set_yaxis_logscale();
      g.set_xaxis_label(combine("1/%s\\S%d",name,(int)pow2));
      g.set_line_style(line_style_t::NO_LINE);
      g.set_settype(settype_t::XYDY);
      g.set_xaxis_min_max(0.001, 1.1);
    }
  
  auto col_scheme=[]() -> vector<grace::color_t>
    {
      using namespace grace;
      return {RED,BLUE,VIOLET,GREEN4,CYAN,MAGENTA,ORANGE,INDIGO,MAROON,TURQUOISE};
    }();
      
  ///loop on different slices
  vector<size_t> comps(3);
  comps[0]=idiag;
  for(size_t jdiv=dim1_min;jdiv<dim1.size();jdiv++)
    {
      comps[idim1]=jdiv;
      double xmin=1e300,xmax=-1e300;
      
      for(auto& p : {&plot_ave,&plot_err})
	{
	  p->new_data_set();
      	  p->set_all_colors(col_scheme[jdiv%col_scheme.size()]);
	}
      
      for(size_t hdiv=dim2_min;hdiv<dim2.size();hdiv++)
	{
	  ///creating ave_err_t variables
	  comps[idim2]=hdiv;
	  const int i=ind(comps);
	  const ave_err_t slop=range_ave_stddev(slopes[i]);
	  const ave_err_t er=range_ave_stddev(errors[i]);
	  
	  ///writing grace files
	  const double x=1.0/pow(dim2[hdiv],pow2);
	  plot_ave.write_ave_err(x,slop);
	  plot_err.write_ave_err(x,er);
	  
	  xmin=min(xmin,x);
	  xmax=max(xmax,x);
	}
      
      plot_err.write_line(bind(fun,std::placeholders::_1,1.0/pow(dim1[jdiv],pow1),pow2),xmin/2.0,xmax*2.0,1001);
      plot_err.set_all_colors(col_scheme[jdiv%col_scheme.size()]);
      
      for(auto& p : {&plot_ave,&plot_err})
	p->set_legend(combine("%s=%d",ind.name(idim1).c_str(),dim1[jdiv]));
    }
}

int main(int narg,char **arg)
{
  //if argument was given, the input is read from pars
  string name="input.txt";
  if(narg>=2) name=arg[1];
  raw_file_t input(name,"r");
  
  cout.precision(16);
  
  T=input.read<int>("T");
  TH=T/2;
  L=input.read<int>("L");
  range_t conf_range;
  input.expect("ConfRange");
  conf_range.start=input.read<size_t>();
  conf_range.each=input.read<size_t>();
  conf_range.end=input.read<size_t>();
  
  set_njacks(input.read<int>("NJacks"));
  const int nhits=input.read<int>("NHits");
  
  const vector<size_t> confs=determine_list_of_confs(conf_range);
  
  //fix on the basis of njacks
  const int nposs_confs=confs.size();
  cout<<"Npossible confs: "<<nposs_confs<<endl;
  const int clust_size=nposs_confs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  vector <int> div_clust_size=divisors(clust_size);
  
  //computes the list of possible number of confs to be used
  vector<int> nconfs;
  for(size_t i=0;i<div_clust_size.size();i++)
    nconfs.push_back(div_clust_size[i]*njacks);
  cout<<"Nconfs: "<<nconfs.back()<<endl;
  
  //computes the list of possible number of hits to be used
  cout<<"Nhits "<<nhits<<endl;
  vector <int> div_nhits=divisors(nhits);
  
  ///checking accessibility of data file. If unaccessible, acquisition routine will start in order to generate it.
  const char data_path[]="plots/data.dat";
  if(not file_exists(data_path))
    {
      cout<<"Data set, absent or unaccessible. Starting generation routine"<<endl;
      
      ///opening output file
      raw_file_t data(data_path,"w");
      
      ///creating indexed dataset for EU5 and EU6
      index_t ind56({{"conf",nconfs.back()},{"diag",2},{"row",(nhits*(nhits-1))/2}});
      vector<dcompl_t> EU56_repository;
      for(int iconf=0;iconf<nconfs.back();iconf++)
	for(int idiag=3;idiag<ndiag;idiag++)
	  {
	    vector<dcompl_t> v=read_vector(combine("out/%04d/EU%d_stoch",confs[iconf],diag[idiag]),(nhits*(nhits-1))/2);
	    EU56_repository.insert(EU56_repository.end(),v.begin(),v.end());
	  }
      ///cycles over configurations numbers, configurations ranges, hits numbers and hits ranges
      for(size_t jdiv=0;jdiv<nconfs.size();jdiv++)
	for(int jrange=0;jrange<nconfs.back()/nconfs[jdiv];jrange++)
	  for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++)
	    {
	      // We can't calculate EU5 and 6 with one single hit, so we avoid to do useless fit
	      const size_t diagmax=(hdiv==0)?3:ndiag;
	      
	      for(int hrange=0;hrange<div_nhits.back()/div_nhits[hdiv];hrange++)
		{
		  /// Allocate diagrams and initializes it
		  array<array<djack_t,2>,ndiag> EU;
		  djvec_t pion(T);
		  array<djvec_t,ndiag> pion_EU;
		  pion_EU.fill(djvec_t{(size_t)T});
		  for(auto &E : EU)
		    for(auto &Eri : E)
		      Eri=0.0;
		  
		  for(size_t iconf=0;iconf<(size_t)nconfs[jdiv];iconf++)
		    {
		      /// Finds the jack index
		      int ijack=iconf/div_clust_size[jdiv];
		      /// Computes the configuration id
		      const int conf=confs[iconf+(njacks*jrange)];
		      
		      /// Reads all diagrams
		      array<vector<dcompl_t>,ndiag> EU_stoch;
		      vector<dcompl_t> pion_stoch=read_vector(combine("out/%04d/mes_contr_Pion",conf),T,5);
		      array<dcompl_t,ndiag> EU_tot;
		      for(int t=0;t<T;t++)
			pion[t][ijack]+=pion_stoch[t].real();
		      
		      for(size_t idiag=0;idiag<(size_t)diagmax;idiag++)
			{
			  if (idiag<3)
			    EU_stoch[idiag]=read_vector(combine("out/%04d/EU%d_stoch",conf,diag[idiag]),div_nhits[hdiv],hrange*div_nhits[hdiv]);
			  ///EU5 and EU6 needs to be read apart because of their structure
			  if(idiag>2)
			    {
			      //routine to extract data with the right pattern
			      size_t count=0;
			      for(int i=hrange*(div_nhits[hdiv])+1;i<(hrange+1)*(div_nhits[hdiv]);i++)
				{
				  count++;
				  for(size_t j=(i*(i+1))/2-count;j<(size_t)((i*(i+1))/2);j++)
				    EU_stoch[idiag].emplace_back(EU56_repository[ind56({iconf+(njacks*jrange),idiag-3,j})]);
				}
			    }
			  // Computes the average over stochastich estimates
			  EU_tot[idiag]=average(EU_stoch[idiag]);
			  
			  // Copy real and imaginary part of the disconnected diagram
			  for(int ri=0;ri<2;ri++)
			    EU[idiag][ri][ijack]+=((double*)&(EU_tot[idiag]))[ri];
			  
			  // Computes the connected and disconnected X connected
			  for(int t=0;t<T;t++)
			    pion_EU[idiag][t][ijack]+=pion_stoch[t].real()*EU_tot[idiag].real();
			}
		    }
		  
		  // Prepares the jackknives of disconnected diagrams and print it
		  for(size_t idiag=0;idiag<diagmax;idiag++)
		    for(int ri=0;ri<2;ri++)
		      EU[idiag][ri].clusterize(div_clust_size[jdiv]);
		  
		  pion.clusterize(div_clust_size[jdiv]).symmetrize();
		  pion.ave_err().write(combine("plots/%03d/pion.xmg",nconfs[jdiv]));
		  djvec_t SL(diagmax);
		  
		  //Creates the jackknife of the connected and disconnected diagrams and plots them
		  for(size_t idiag=0;idiag<diagmax;idiag++)
		    {
		      pion_EU[idiag].clusterize(div_clust_size[jdiv]).symmetrize();
		      pion_EU[idiag].ave_err().write(combine("plots/%03d/%03d/pion_EU%d.xmg",nconfs[jdiv],div_nhits[hdiv],diag[idiag]));
		      
		      /// Subtract the fully disconnected part and plots it
		      djvec_t pion_EU_sub=pion_EU[idiag]-pion*EU[idiag][0];
		      pion_EU_sub.ave_err().write(combine("plots/%03d/%03d/EU%d_rangeconfs_%03d_rangehits_%03d_sub.xmg",nconfs[jdiv],div_nhits[hdiv],diag[idiag],jrange,hrange));
		      
		      /// Computes the ratio with the purely connected and plots it
		      djvec_t pion_EU_rat=pion_EU_sub/pion;
		      
		      djack_t Z2,M,DZ2_fr_Z2;
		      two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL[idiag],pion,pion_EU_sub,TH,12,TH,combine("plots/%03d/pion.xmg",nconfs[jdiv]),combine("plots/%03d/%03d/EU%d_rangeconfs_%03d_rangehits_%03d_rat.xmg",nconfs[jdiv],div_nhits[hdiv],diag[idiag],jrange,hrange));
		    }
		  cout<<endl;
		  cout << "confs: " << nconfs[jdiv] << "; range confs: " << jrange << "; hits: " << div_nhits[hdiv] << "; range hits: " << hrange <<endl;
		  
		  ///writing on data file
		  for(size_t idiag=0;idiag<diagmax;idiag++)
		    {
		      data.bin_write(SL[idiag].err());
		      data.bin_write(SL[idiag].ave());
		      cout<<"EU"<<diag[idiag]<<": "<<SL[idiag].ave_err()<<endl;
		    }
		}
	    }
      cout<< "dataset generated, starting analysis"<<endl;
    }
  else
    cout<<"dataset checked, starting analysis"<<endl;
  
  ///opening data file
  raw_file_t data(data_path,"r");
  
  ///performing averages and errors and putting everything in output files.
  index_t ind({{"diagrams",ndiag},{"nhits",div_nhits.size()},{"nconfs",nconfs.size()}});
  vector<vector<double>> slopes(ind.max(),vector<double>(0)),errors(ind.max(),vector<double>(0));
  for(size_t jdiv=0;jdiv<nconfs.size();jdiv++)
    for(int jrange=0;jrange<nconfs.back()/nconfs[jdiv];jrange++)
      for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++)
	{
	  const size_t diagmax=(hdiv==0)?3:ndiag;
	  
	  for(int hrange=0;hrange<div_nhits.back()/div_nhits[hdiv];hrange++)
	    for(size_t idiag=0;idiag<diagmax;idiag++)
	      for(auto& es : {&errors,&slopes})
		(*es)[ind({idiag,hdiv,jdiv})].push_back(data.bin_read<double>());
	}
  
  ///loop on different slices
  const int pow_nh[5]={1,1,1,2,2};
  for(size_t idiag=0;idiag<ndiag;idiag++)
    {
      int seed=0;
      
      cout<<"DIAG "<<diag[idiag]<<endl;
      int idiv_nhits_min=0;
      if(idiag>2)
	idiv_nhits_min=1;
      
      /////////////////////////////////////////////////////////////////
      
      plan_fit_data_t<djack_t> plan_fit_data;
      const int nx=2;
      for(size_t idiv_nhits=idiv_nhits_min;idiv_nhits<div_nhits.size();idiv_nhits++)
	for(size_t idiv_nconfs=0;idiv_nconfs<nconfs.size();idiv_nconfs++)
	  {
	    const int i=ind({idiag,idiv_nhits,idiv_nconfs});
	    
	    vector<double> x(nx);
	    x[0]=1.0/(nconfs[idiv_nconfs]*Ndd(div_nhits[idiv_nhits],idiag));
	    x[1]=1.0/nconfs[idiv_nconfs];
	    
	    const ave_err_t ae=range_ave_stddev(errors[i]);
	    if(ae.err()>1e-10)
	      {
		djack_t j;
		j.fill_gauss(sqr(ae.ave()),ae.ave()*ae.err()*2,seed++);
		plan_fit_data.push_back(make_tuple(x,j));
	      }
	  }
      
      //! results for the fit
      const djvec_t res_plan=plan_fit(plan_fit_data);
      cout<<"res_plan of diagram EU"<<diag[idiag]<<": "<<res_plan.ave_err()<<endl;
      const djvec_t temp=sqrt(abs(res_plan));
      /////////////////////////////////////////////////////////////////
      
      djvec_t res(2);
      jack_fit_t jack_fit;
      // const size_t ipconfs=jack_fit.add_fit_par_limits(res[0],"pconfs",21,0.1,0.0,40);
      // const size_t iphits=jack_fit.add_fit_par_limits(res[1],"phits",11,0.1,0.0,40);
      const size_t ipconfs=(idiag!=3)?jack_fit.add_fit_par(res[0],"pconfs",temp[0].ave(),temp[0].err()):jack_fit.add_fit_par_limits(res[0],"pconfs",0.1,0.1,0.0,1e3);
      const size_t iphits=(idiag!=3)?jack_fit.add_fit_par(res[1],"phits",temp[1].ave(),temp[1].err()):jack_fit.add_fit_par_limits(res[1],"phits",0.1,0.1,0.0,1e3);
      
      for(size_t idiv_nhits=idiv_nhits_min;idiv_nhits<div_nhits.size();idiv_nhits++)
	for(size_t idiv_nconfs=0;idiv_nconfs<nconfs.size();idiv_nconfs++)
	  {
	    const int i=ind({idiag,idiv_nhits,idiv_nconfs});
	    
	    // vector<double> x(nx);
	    // x[0]=1.0/(nconfs[idiv_nconfs]*pow(div_nhits[idiv_nhits],pow_nh[idiag]));
	    // x[1]=1.0/nconfs[idiv_nconfs];
	    
	    const ave_err_t ae=range_ave_stddev(errors[i]);
	    if(ae.err()>1e-6*fabs(ae.ave()))
	      {
		djack_t j;
		j.fill_gauss(sqr(ae.ave()),ae.ave()*ae.err()*2,seed++);
		//j.fill_gauss((20/nconfs[idiv_nconfs]+0.1/div_nhits[idiv_nconfs]),0.00001,seed++);
		//plan_fit_data.push_back(make_tuple(x,j));
		
		jack_fit.add_point(//numerical data
				   [j// ,i,ae,idiv_nhits,idiv_nconfs
				    ]
				   (const vector<double> &p,int iel)
				   {
				     // static set<pair<int,int>> l;
				     // if(l.find({iel,i})==l.end())
				     //   {
				     // 	 l.insert({iel,i});
				     // 	 cout<<iel<<" "<<i<<" "<<j[iel]<<"     "<<ae<<"    idiv_nhit: "<<idiv_nhits<<"    idiv_nconfs: "<<idiv_nconfs<<endl;
				     //   }
				     return j[iel];
				   },
				   //ansatz
				   [ipconfs,nconfs,idiv_nconfs,div_nhits,idiv_nhits,idiag,iphits]
				   (const vector<double> &p,int iel)
				   {
				     const double A=sqr(p[ipconfs]);
				     const double B=sqr(p[iphits]);
				     return
				       A/(nconfs[idiv_nconfs]*Ndd(div_nhits[idiv_nhits],idiag))+
				       B/nconfs[idiv_nconfs];
				   },
				   
				   //for covariance/error
				   j.err());
	      }
	  }
      
      //! results for the fit
      // const djvec_t res=plan_fit(plan_fit_data);
      jack_fit.fit();
      
      
      //! fit ansatz
      auto f=[&res](const double inv_nconfs,const double x,const int pow)
	{
	  const djack_t A=sqr(res[0]);
	  const djack_t B=sqr(res[1]);
	  const double inv_Ndd=(pow==1)?x:(x/(1-sqrt(x)));
	  //if(pow!=1) cout<<pow<<" "<<x<<" "<<inv_Ndd<<endl;
	  return djack_t(sqrt((A*inv_Ndd+B)*inv_nconfs)).ave();
	};
      
      //! same ansatz with swapped arguments
      auto f_swapped=bind(f,placeholders::_2,placeholders::_1,placeholders::_3);
      
      cout<<res.ave_err()<<endl;
      const djack_t r=res[0]/res[1];
      const djack_t lim=pow(sqr(r),1.0/pow_nh[idiag]);
      cout<<"the 'nhits'-associated error is no longer dominating for nhits which equals "<<smart_print(lim.ave_err())<<endl;
      
      slice_plot(div_nhits,idiv_nhits_min,nconfs,0,idiag,slopes,errors,ind,1,2,f,pow_nh[idiag],1);
      slice_plot(nconfs,0,div_nhits,idiv_nhits_min,idiag,slopes,errors,ind,2,1,f_swapped,1,pow_nh[idiag]);
    }
  
  return 0;
}
