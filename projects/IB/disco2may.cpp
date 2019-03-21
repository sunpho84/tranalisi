#include <tranalisi.hpp>

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
void slice_plot(const vector<int>& dim1,const vector<int>& dim2,const size_t idiag,const vector<vector<double>>& slopes,const vector<vector<double>>& errors,const index_t& ind,int idim1,int idim2,const F& fun)
{
  const char* name=ind.name(idim2).c_str();
  grace_file_t plot_ave(combine("plots/fit_results/EU%d_average_vs_%s",diag[idiag],name));
  grace_file_t plot_err(combine("plots/fit_results/EU%d_errors_vs_%s",diag[idiag],name));
  bool ytitle[5]={1,1,0,0,0};
  
  if(ytitle[idiag]==1)
    {
      plot_ave.set_yaxis_label("averages");
      plot_err.set_yaxis_label("errors");
    }
  else
    {
      plot_ave.set_yaxis_label("averages (a*M)");
      plot_err.set_yaxis_label("errors (a*M)");
    }
  for(grace_file_t* _g : {&plot_ave,&plot_err})
    {
      grace_file_t& g=(*_g);
      
      using namespace grace;
      
      g.set_color_scheme({ORANGE,RED,BLUE,GREEN4});
      
      g.set_title(combine("Errors vs %s",name));
      g.set_subtitle(combine("EU%d",diag[idiag]));
      g.set_xaxis_logscale();
      g.set_yaxis_logscale();
      g.set_xaxis_label(combine("1/%s",name));
      g.set_line_style(line_style_t::NO_LINE);
      g.set_settype(settype_t::XYDY);
      g.set_xaxis_min_max(0.001, 1.1);
    }
  
  ///loop on different slices
  vector<size_t> comps(3);
  comps[0]=idiag;
  for(size_t jdiv=0;jdiv<dim1.size();jdiv++)
    {
      comps[idim1]=jdiv;
      double xmin=1e300,xmax=-1e300;
      for(size_t hdiv=0;hdiv<dim2.size();hdiv++)
	{
	  ///creating ave_err_t variables
	  comps[idim2]=hdiv;
	  const int i=ind(comps);
	  const ave_err_t slop=range_ave_stddev(slopes[i]);
	  const ave_err_t er=range_ave_stddev(errors[i]);
	  
	  ///writing grace files
	  const double x=1.0/dim2[hdiv];
	  plot_ave.write_ave_err(x,slop);
	  plot_err.write_ave_err(x,er);
	  
	  xmin=min(xmin,x);
	  xmax=max(xmax,x);
	}
      
      plot_err.write_line(bind(fun,std::placeholders::_1,1.0/dim1[jdiv]),xmin/10.0,xmax*10.0);
      
      for(auto& _p : {&plot_ave,&plot_err})
	{
	  grace_file_t& p=(*_p);
	  p.set_legend(combine("%s=%d",ind.name(idim1).c_str(),dim1[jdiv]));
	  p.new_data_set();
	}
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
      
      ///cycles over configurations numbers, configurations ranges, hits numbers and hits ranges
      for(size_t jdiv=0;jdiv<nconfs.size();jdiv++)
	for(int jrange=0;jrange<nconfs.back()/nconfs[jdiv];jrange++)
	  for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++)
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
		
		for(int iconf=0;iconf<nconfs[jdiv];iconf++)
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
		    
		    for(int idiag=0;idiag<ndiag;idiag++)
		      {
			EU_stoch[idiag]=read_vector(combine("out/%04d/EU%d_stoch",conf,diag[idiag]),div_nhits[hdiv],hrange*div_nhits[hdiv]);
			
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
		for(int idiag=0;idiag<ndiag;idiag++)
		  for(int ri=0;ri<2;ri++)
		    EU[idiag][ri].clusterize(div_clust_size[jdiv]);
		
		pion.clusterize(div_clust_size[jdiv]).symmetrize();
		pion.ave_err().write(combine("plots/%03d/pion.xmg",nconfs[jdiv]));
		djvec_t SL(ndiag);
		
		// Creates the jackknife of the connected and disconnected diagrams and plots them
		for(int idiag=0;idiag<ndiag;idiag++)
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
		for(size_t idiag=0;idiag<ndiag;idiag++)
		  {
		    data.bin_write(SL[idiag].err());
		    data.bin_write(SL[idiag].ave());
		    cout<<"EU"<<diag[idiag]<<": "<<SL[idiag].ave_err()<<endl;
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
	for(int hrange=0;hrange<div_nhits.back()/div_nhits[hdiv];hrange++)
	  for(size_t idiag=0;idiag<ndiag;idiag++)
	    for(auto& es : {&errors,&slopes})
	      (*es)[ind({idiag,hdiv,jdiv})].push_back(data.bin_read<double>());
  
  ///loop on different slices
  for(size_t idiag=0;idiag<ndiag;idiag++)
    {
      const int nx=3;
      plan_fit_data_t<djack_t> plan_fit_data;
      int seed=0;
      
      cout<<"DIAG "<<diag[idiag]<<endl;
      
      for(size_t idiv_nhits=0;idiv_nhits<div_nhits.size();idiv_nhits++)
	for(size_t idiv_nconfs=0;idiv_nconfs<nconfs.size();idiv_nconfs++)
	  {
	    ///creating ave_err_t variables
	    const int i=ind({idiag,idiv_nhits,idiv_nconfs});
	    
	    vector<double> x(nx);
	    x[0]=1.0;
	    x[1]=1.0/(nconfs[idiv_nconfs]*div_nhits[idiv_nhits]);
	    x[2]=1.0/nconfs[idiv_nconfs];
	    
	    const ave_err_t ae=range_ave_stddev(errors[i]);
	     if(ae.err()>1e-10)
	      {
		djack_t j;
		j.fill_gauss(sqr(ae.ave()),ae.ave()*ae.err()*2,seed++);
		//j.fill_gauss((20/nconfs[idiv_nconfs]+0.1/div_nhits[idiv_nconfs]),0.00001,seed++);
		plan_fit_data.push_back(make_tuple(x,j));
	      }
	  }
      
      const djvec_t res=plan_fit(plan_fit_data,true);
      
      //! fit ansatz
      auto f=[&res](double inv_nconfs,double inv_nhits)
	{
	  return djack_t(sqrt(res[0]+(res[1]*inv_nhits+res[2])*inv_nconfs)).ave();
	};
      
      //! same ansatz with swapped arguments
      auto f_swapped=bind(f,placeholders::_2,placeholders::_1);
      
      //cout<<res.ave_err()<<endl;
      
      slice_plot(div_nhits,nconfs,idiag,slopes,errors,ind,1,2,f);
      slice_plot(nconfs,div_nhits,idiag,slopes,errors,ind,2,1,f_swapped);
    }
  
  return 0;
}
