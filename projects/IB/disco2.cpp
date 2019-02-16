#include <tranalisi.hpp>

int T,L;
int TH;

vector<dcomplex> read_vector(const string &path,int n,int nskip=0)
{
  raw_file_t data(path,"r");
  
  cout<<"Reading from "<<path<<" "<<n<<" lines after skipping "<<nskip<<" words"<<endl;
  for(int iskip=0;iskip<nskip;iskip++)
    data.skip_line();
  
  vector<dcomplex> d;
  
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
  return accumulate(v.begin(),v.end(),dcomplex{0,0})/(double)v.size();
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
  
  vector<size_t> confs_p=get_existing_paths_in_range("out/%04d/mes_contr_Pion",conf_range);
  vector<size_t> confs_1=get_existing_paths_in_range("out/%04d/EU1_stoch",conf_range);
  vector<size_t> confs;
  for(auto &_p : confs_p)
    for(auto & _1 : confs_1)
      if(_p==_1)
	{
	  cout<<_p<<endl;
	  confs.push_back(_p);
	}
  
  set_njacks(input.read<int>("NJacks"));
  
  const int nposs_confs=confs.size();
  cout<<"Npossible confs: "<<nposs_confs<<endl;
  const int clust_size=nposs_confs/njacks;
  cout<<"Cluster size: "<<clust_size<<endl;
  const int nconfs=clust_size*njacks;
  cout<<"Nconfs: "<<nconfs<<endl;
  
  const int nhits=input.read<int>("NHits");
  
  /// Allocate diagrams and initializes it
  djack_t EU1[2]={0.0,0.0};
  djack_t EU2[2]={0.0,0.0};
  djack_t EU4[2]={0.0,0.0};
  djack_t EU5[2]={0.0,0.0};
  djack_t EU6[2]={0.0,0.0};
  djvec_t pion(T);
  djvec_t pion_EU1(T);
  djvec_t pion_EU2(T);
  djvec_t pion_EU4(T);
  djvec_t pion_EU5(T);
  djvec_t pion_EU6(T);
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      /// Finds the jack index
      int ijack=iconf/clust_size;
      
      /// Computes the configuration id
      const int conf=confs[iconf];
      
      /// Reads all diagrams
      vector<dcomplex> EU1_stoch=read_vector(combine("out/%04d/EU1_stoch",conf),nhits);
      vector<dcomplex> EU2_stoch=read_vector(combine("out/%04d/EU2_stoch",conf),nhits);
      vector<dcomplex> EU4_stoch=read_vector(combine("out/%04d/EU4_stoch",conf),nhits);
      vector<dcomplex> EU5_stoch=read_vector(combine("out/%04d/EU5_stoch",conf),nhits*(nhits-1)/2);
      vector<dcomplex> EU6_stoch=read_vector(combine("out/%04d/EU6_stoch",conf),nhits*(nhits-1)/2);
      vector<dcomplex> pion_stoch=read_vector(combine("out/%04d/mes_contr_Pion",conf),T,5);
      
      /// Computes the average over stochastich estimates
      dcomplex EU1_tot=average(EU1_stoch);
      dcomplex EU2_tot=average(EU2_stoch);
      dcomplex EU4_tot=average(EU4_stoch);
      dcomplex EU5_tot=average(EU5_stoch);
      dcomplex EU6_tot=average(EU6_stoch);
      
      /// Copy real and imaginary part of the disconnected diagram
      for(int ri=0;ri<2;ri++)
	{
	  EU1[ri][ijack]+=((double*)&EU1_tot)[ri];
	  EU2[ri][ijack]+=((double*)&EU2_tot)[ri];
	  EU4[ri][ijack]+=((double*)&EU4_tot)[ri];
	  EU5[ri][ijack]+=((double*)&EU5_tot)[ri];
	  EU6[ri][ijack]+=((double*)&EU6_tot)[ri];
	}
      
      /// Computes the connected and disconnected X connected
      for(int t=0;t<T;t++)
	{
	  pion[t][ijack]+=pion_stoch[t].real();
          pion_EU1[t][ijack]+=pion_stoch[t].real()*EU1_tot.real();
          pion_EU2[t][ijack]+=pion_stoch[t].real()*EU2_tot.real();
          pion_EU4[t][ijack]+=pion_stoch[t].real()*EU4_tot.real();
          pion_EU5[t][ijack]+=pion_stoch[t].real()*EU5_tot.real();
          pion_EU6[t][ijack]+=pion_stoch[t].real()*EU6_tot.real();
	}
    }
  
  /// Prepares the jackknives of disconnected diagrams and print it
  for(int ri=0;ri<2;ri++)
    {
      EU1[ri].clusterize(clust_size);
      EU2[ri].clusterize(clust_size);
      EU4[ri].clusterize(clust_size);
      EU5[ri].clusterize(clust_size);
      EU6[ri].clusterize(clust_size);
      
      cout<<"EU1: "<<EU1[ri].ave_err()<<endl;
      cout<<"EU2: "<<EU2[ri].ave_err()<<endl;
      cout<<"EU4: "<<EU4[ri].ave_err()<<endl;
      cout<<"EU5: "<<EU5[ri].ave_err()<<endl;
      cout<<"EU6: "<<EU6[ri].ave_err()<<endl;
    }
  
  /// Creates the jackknife of the connected and disconnected diagrams and plots them
  pion.clusterize(clust_size).symmetrize();
  pion_EU1.clusterize(clust_size).symmetrize();
  pion_EU2.clusterize(clust_size).symmetrize();
  pion_EU4.clusterize(clust_size).symmetrize();
  pion_EU5.clusterize(clust_size).symmetrize();
  pion_EU6.clusterize(clust_size).symmetrize();
  
  pion.ave_err().write("plots/pion.xmg");
  pion_EU1.ave_err().write("plots/pion_EU1.xmg");
  pion_EU2.ave_err().write("plots/pion_EU2.xmg");
  pion_EU4.ave_err().write("plots/pion_EU4.xmg");
  pion_EU5.ave_err().write("plots/pion_EU5.xmg");
  pion_EU6.ave_err().write("plots/pion_EU6.xmg");
  
  /// Subtract the fully disconnected part and plots it
  djvec_t pion_EU1_sub=pion_EU1-pion*EU1[0];
  djvec_t pion_EU2_sub=pion_EU2-pion*EU2[0];
  djvec_t pion_EU4_sub=pion_EU4-pion*EU4[0];
  djvec_t pion_EU5_sub=pion_EU5-pion*EU5[0];
  djvec_t pion_EU6_sub=pion_EU6-pion*EU6[0];
  pion_EU1_sub.ave_err().write("plots/EU1_sub.xmg");
  pion_EU2_sub.ave_err().write("plots/EU2_sub.xmg");
  pion_EU4_sub.ave_err().write("plots/EU4_sub.xmg");
  pion_EU5_sub.ave_err().write("plots/EU5_sub.xmg");
  pion_EU6_sub.ave_err().write("plots/EU6_sub.xmg");
  
  /// Computes the ratio with the purely connected and plots it
  djvec_t pion_EU1_rat=pion_EU1_sub/pion;
  djvec_t pion_EU2_rat=pion_EU2_sub/pion;
  djvec_t pion_EU4_rat=pion_EU4_sub/pion;
  djvec_t pion_EU5_rat=pion_EU5_sub/pion;
  djvec_t pion_EU6_rat=pion_EU6_sub/pion;
  
  djack_t Z2,M,DZ2_fr_Z2;
  djack_t SL_1,SL_2,SL_4,SL_5,SL_6;
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL_1,pion,pion_EU1_sub,TH,12,TH,"plots/pion.xmg","plots/EU1_rat.xmg");
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL_2,pion,pion_EU2_sub,TH,12,TH,"plots/pion.xmg","plots/EU2_rat.xmg");
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL_4,pion,pion_EU4_sub,TH,12,TH,"plots/pion.xmg","plots/EU4_rat.xmg");
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL_5,pion,pion_EU5_sub,TH,12,TH,"plots/pion.xmg","plots/EU5_rat.xmg");
  two_pts_with_ins_ratio_fit(Z2,M,DZ2_fr_Z2,SL_6,pion,pion_EU6_sub,TH,12,TH,"plots/pion.xmg","plots/EU6_rat.xmg");
  
  cout<<SL_1.ave_err()<<endl;
  cout<<SL_2.ave_err()<<endl;
  cout<<SL_4.ave_err()<<endl;
  cout<<SL_5.ave_err()<<endl;
  cout<<SL_6.ave_err()<<endl;
  
  // pion_EU1_rat.ave_err().write("plots/EU1_rat.xmg");
  // pion_EU2_rat.ave_err().write("plots/EU2_rat.xmg");
  // pion_EU4_rat.ave_err().write("plots/EU4_rat.xmg");
  // pion_EU5_rat.ave_err().write("plots/EU5_rat.xmg");
  // pion_EU6_rat.ave_err().write("plots/EU6_rat.xmg");
  
  return 0;
}
