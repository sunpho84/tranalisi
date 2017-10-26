#include <tranalisi.hpp>

#include <set>

//master_fprintf(fout,"%d\t",iconf);
//master_fprintf(fout,"%d\t",icopy);
using tmeas_t=size_t;
using flav_t=size_t;
using op_t=size_t;
			// master_fprintf(fout,"%+16.16lg\t",plaq);
			// master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tot_charge,tot_charge2);
			// master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",spinpol[RE],spinpol[IM]);
			// master_fprintf(fout,"%+16.16lg" "\t" "%+16.16lg\t",tens[RE],tens[IM]);

using obs_t=valarray<double>;
enum{PLAQ,TOT_CHARGE,TOT_CHARGE2,SPI_RE,SPI_IM,TENS_RE,TENS_IM};
const size_t OBS_SIZE=7;

//! returns the position of the tag in the list
size_t get_op_pos(const set<op_t> &op_list,const op_t &tag)
{
  auto op_it=op_list.find(tag);
  if(op_it==op_list.end()) CRASH("unable to find operator %zu",tag);
  return distance(op_list.begin(),op_it);
}

size_t ntimes;
size_t ncopies;
size_t nflavs;
size_t nops;
size_t clust_size;
size_t nconfs;

//! filter an observable
djack_t get_obs(const vector<obs_t> &obs,const index_t &ind,const size_t itime,const size_t iflav,const size_t iop,size_t iobs)
{
  djack_t out;
  out=0.0;
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    for(size_t icopy=0;icopy<ncopies;icopy++)
      {
	size_t iclust=iconf/clust_size;
	out[iclust]+=obs[ind({iconf,icopy,itime,iflav,iop})][iobs];
      }
  
  out.clusterize();
  
  return out;
}

int main()
{
  set_njacks(3);
  
  set<tmeas_t> tmeas_list;
  set<flav_t> flav_list;
  set<op_t> op_list;
  
  char pollo_path[]="pollo";
  ifstream pollo(pollo_path);
  if(not pollo.good()) CRASH("opening %s",pollo_path);
  
  ncopies=0;
  
  size_t conf,copy,tmeas,flav;
  string op_str;
  obs_t obs(OBS_SIZE);
  vector<obs_t> obs_list;
  while(pollo>>conf>>copy>>tmeas>>flav>>op_str>>obs[PLAQ]>>obs[TOT_CHARGE]>>obs[TOT_CHARGE2]>>obs[SPI_RE]>>obs[SPI_IM]>>obs[TENS_RE]>>obs[TENS_IM])
    {
      size_t op;
      const char format_str[]="%zu,0";
      if(sscanf(op_str.c_str(),format_str,&op)!=1) CRASH("parsing string %s, expecting something in the format",op_str.c_str(),format_str);
      
      tmeas_list.insert(tmeas);
      flav_list.insert(flav);
      op_list.insert(op);
      
      ncopies=max(ncopies,copy+1);
      obs_list.push_back(obs);
    }
  
  ntimes=tmeas_list.size();
  nflavs=flav_list.size();
  nops=op_list.size();
  size_t nconfs_max=obs_list.size()/(ntimes*nflavs*nops);
  clust_size=nconfs_max/njacks;
  nconfs=clust_size*njacks;
  
  cout<<"NTimes: "<<ntimes<<endl;
  cout<<"NFlavs: "<<nflavs<<endl;
  cout<<"Nops: "<<nops<<endl;
  cout<<"NCopies: "<<ncopies<<endl;
  cout<<"NConfs_max: "<<nconfs_max<<endl;
  cout<<"Clust_size: "<<clust_size<<endl;
  cout<<"Nconfs: "<<nconfs<<endl;
  
  const size_t sigma_xy_tag=6,i_sigma_xy=get_op_pos(op_list,sigma_xy_tag);
  const size_t sigma_zt_tag=9,i_sigma_zt=get_op_pos(op_list,sigma_zt_tag);
  
  const index_t ind({{"conf",nconfs},{"copy",ncopies},{"tmeas",ntimes},{"flav",nflavs},{"op",nops}});
  
  cout<<"Operator Sigma_xy pos: "<<i_sigma_xy<<endl;
  cout<<"Operator Sigma_zt pos: "<<i_sigma_zt<<endl;

  enum{UP,DW,ST};
  for(size_t iflav=0;iflav<nflavs;iflav++)
    {
      cout<<"-------------------------------------------------------"<<endl;
      cout<<"iflav: "<<iflav<<endl;
      djack_t sigma_xy=get_obs(obs_list,ind,ntimes-1,iflav,i_sigma_xy,TENS_IM);
      djack_t Qsigma_zt=get_obs(obs_list,ind,ntimes-1,iflav,i_sigma_zt,SPI_IM);
      djack_t rat=Qsigma_zt/sigma_xy*1000;
      cout<<"rat: "<<smart_print(rat)<<endl;
      djack_t Q=get_obs(obs_list,ind,ntimes-1,iflav,i_sigma_xy,TOT_CHARGE);
      cout<<"Q: "<<smart_print(Q)<<endl;
      size_t glb_vol=110592;
      djack_t Q2=get_obs(obs_list,ind,ntimes-1,iflav,i_sigma_xy,TOT_CHARGE2)/glb_vol;
      cout<<"Q2: "<<smart_print(Q2)<<endl;
      djack_t C_f=Qsigma_zt/(sigma_xy*sqrt(Q2));
      cout<<"C_f: "<<smart_print(C_f)<<endl;
    }
  
  return 0;
}
