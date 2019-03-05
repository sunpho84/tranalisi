#include <tranalisi.hpp>

int T,L;
int TH;

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
vector <int> divisors(const int n)
{
  vector <int> b;
  for(int idiv=1;idiv<n+1;idiv++) if (n%idiv==0) b.push_back(idiv);
  return b;
}
///checking input files accessibility
bool checkmethod(string str)
{
  ifstream f(str);
  return f.good();
}

int main(int narg,char **arg)
{
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
  
  vector<size_t> confs_p=get_existing_paths_in_range("out/%04d/mes_contr_Pion",conf_range,SILENT);
  vector<size_t> confs_1=get_existing_paths_in_range("out/%04d/EU1_stoch",conf_range,SILENT);
  vector<size_t> confs;
  for(auto &_p : confs_p)
    for(auto & _1 : confs_1)
      if(_p==_1)
	{
	  //	  cout<<_p<<endl;
	  confs.push_back(_p);
	}
  
  set_njacks(input.read<int>("NJacks"));
  
  const int nposs_confs=confs.size();
  cout<<"Npossible confs: "<<nposs_confs<<endl;
  const int clust_size=nposs_confs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  vector <int> div_clust_size=divisors(clust_size);
  
  vector <int> nconfs;
  for(int i=0;i<div_clust_size.size();i++) nconfs.push_back(div_clust_size[i]*njacks);
  cout<<"Nconfs: "<<nconfs.back()<<endl;
  const int nhits=input.read<int>("NHits");
  cout<<"Nhits "<<nhits<<endl;
  vector <int> div_nhits=divisors(nhits);
  const int ndiag=5;
  enum {EU1,EU2,EU4,EU5,EU6};
  const int diag[ndiag]={1,2,4,5,6};
  
  
  ///checking accessibility of data file. If unaccessible, acquisition routine will start in order to generate it.
  if(!checkmethod("plots/data.dat")){
    cout<<"Data set incomplete, absent or unaccessible. Starting generation routine"<<endl;
    ///opening output file
    raw_file_t data("plots/data.dat","w");
    ///cycles over configurations numbers, configurations ranges, hits numbers and hits ranges
    for(size_t jdiv=0;jdiv<nconfs.size();jdiv++){
      for(size_t jrange=0;jrange<nconfs.back()/nconfs[jdiv];jrange++){
	for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++){
	  for(size_t hrange=0;hrange<div_nhits.back()/div_nhits[hdiv];hrange++){
	    
	    
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
		{
		  EU[idiag][ri].clusterize(div_clust_size[jdiv]);
		  
		  //cout<<"EU"<<diag[idiag]<<": "<<EU[idiag][ri].ave_err()<<endl;
		}
	    
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
	}
      }
    }
    cout<< "dataset generated, starting analysis"<<endl;
  }
  else {cout<<"dataset checked, starting analysis"<<endl;}
  ///opening data file
  raw_file_t data("plots/data.dat","r");
  ///opening output
  vector<grace_file_t> plot_conf_ave(ndiag), plot_conf_err(ndiag), plot_hits_ave(ndiag), plot_hits_err(ndiag);
  for(int idiag=0;idiag<ndiag;idiag++)
	{
	  plot_conf_ave[idiag].open(combine("plots/fit_results/EU%d_average_vs_nconfs",diag[idiag]));
	  plot_conf_err[idiag].open(combine("plots/fit_results/EU%d_errors_vs_nconfs",diag[idiag]));
	  plot_hits_ave[idiag].open(combine("plots/fit_results/EU%d_average_vs_nhits",diag[idiag]));
	  plot_hits_err[idiag].open(combine("plots/fit_results/EU%d_errors_vs_nhits",diag[idiag]));
	}
  ///performing averages and errors and putting everything in output files.
  double dat;
  index_t ind({{"diagrams",ndiag},{"nhits",div_nhits.size()}});
  vector<vector<double>> slopes(ind.max(),vector<double>(0)),errors(ind.max(),vector<double>(0));
    for(size_t jdiv=0;jdiv<nconfs.size();jdiv++){
      array<vector<double>,ndiag> slopesconf,errorsconf;
      for(size_t jrange=0;jrange<nconfs.back()/nconfs[jdiv];jrange++){
	for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++){
	  for(size_t hrange=0;hrange<div_nhits.back()/div_nhits[hdiv];hrange++){
	    for(size_t idiag=0;idiag<ndiag;idiag++){
	      data.bin_read(dat);
	      errors[ind({idiag,hdiv})].push_back(dat);
	      errorsconf[idiag].push_back(dat);
	      data.bin_read(dat);
	      slopes[ind({idiag,hdiv})].push_back(dat);
	      slopesconf[idiag].push_back(dat);
	    }
	  }
	}
      }
      ///printing average and errors vs. nconfs
      for(int idiag=0;idiag<ndiag;idiag++)
	{
	  ///generating valarray from vectors
	  valarray<double> m(slopesconf[idiag].data(),slopesconf[idiag].size()),e(errorsconf[idiag].data(),errorsconf[idiag].size());
	  ///creating ave_err_t variables
	  ave_err_t slop, er;
	  slop=range_ave_stddev(m,slopesconf[idiag].size());
	  er=range_ave_stddev(e,errorsconf[idiag].size());
	  ///writing grace files
	  plot_conf_ave[idiag].write_ave_err(nconfs[jdiv],slop);
	  plot_conf_err[idiag].write_ave_err(nconfs[jdiv],er);
	}
    }
    ///printing average and errors vs. 1/nhits
    for(size_t idiag=0;idiag<ndiag;idiag++)
      {
	for(size_t hdiv=0;hdiv<div_nhits.size();hdiv++)
	  {
	    ///generating valarray from vectors
	    valarray<double> m(slopes[ind({idiag,hdiv})].data(),slopes[ind({idiag,hdiv})].size()),e(errors[ind({idiag,hdiv})].data(),errors[ind({idiag,hdiv})].size());
	    ///creating ave_err_t variables
	    ave_err_t slop, er;
	    slop=range_ave_stddev(m,slopes[ind({idiag,hdiv})].size());
	    er=range_ave_stddev(e,errors[ind({idiag,hdiv})].size());
	    ///writing grace files
	    plot_hits_ave[idiag].write_ave_err(1./div_nhits[hdiv],slop);
	    plot_hits_err[idiag].write_ave_err(1./div_nhits[hdiv],er);
	  }
      }
    return 0;
}
