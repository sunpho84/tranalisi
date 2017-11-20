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
size_t nconf_types;

size_t clust_size;
size_t nconfs;

//! filter an observable
djack_t get_obs(const vector<obs_t> &obs,const index_t &ind,const size_t itime,const size_t iflav,const size_t iop,const size_t gauge_ferm_conf,size_t iobs)
{
  djack_t out;
  out=0.0;
  
  for(size_t iconf=0;iconf<nconfs;iconf++)
    for(size_t icopy=0;icopy<ncopies;icopy++)
      {
	size_t iclust=iconf/clust_size;
	size_t ientry=ind({iconf,icopy,itime,iflav,iop,gauge_ferm_conf});
	out[iclust]+=obs[ientry][iobs];
      }
  
  out.clusterize(clust_size);
  out/=ncopies;
  
  return out;
}

int main()
{
  const size_t ext_njacks=15;
  set_njacks(ext_njacks);
  cout<<"NJacks: "<<njacks<<endl;
  
  const size_t L=24,T=8;
  
  set<tmeas_t> tmeas_list;
  set<flav_t> flav_list;
  set<op_t> op_list;
  set<size_t> conf_type_list;
  
  char pollo_path[]="pollo";
  ifstream pollo(pollo_path);
  if(not pollo.good()) CRASH("opening %s",pollo_path);
  
  vector<obs_t> obs_list;
  ncopies=0;
  {
    size_t conf,copy,tmeas,flav;
    string op_str;
    int conf_type;
    obs_t obs(OBS_SIZE);
    while(pollo>>conf>>copy>>tmeas>>flav>>op_str>>conf_type>>obs[PLAQ]>>obs[TOT_CHARGE]>>obs[TOT_CHARGE2]>>obs[SPI_RE]>>obs[SPI_IM]>>obs[TENS_RE]>>obs[TENS_IM])
      {
	size_t op;
	const char format_str[]="%zu,0";
	if(sscanf(op_str.c_str(),format_str,&op)!=1) CRASH("parsing string %s, expecting something in the format",op_str.c_str(),format_str);
	
	tmeas_list.insert(tmeas);
	flav_list.insert(flav);
	op_list.insert(op);
	conf_type_list.insert(conf_type);
	
	ncopies=max(ncopies,copy+1);
	obs_list.push_back(obs);
      }
  }
  
  ntimes=tmeas_list.size();
  nflavs=flav_list.size();
  nops=op_list.size();
  nconf_types=conf_type_list.size();
  size_t nconfs_max=obs_list.size()/(ntimes*nflavs*nops*nconf_types);
  clust_size=nconfs_max/njacks;
  nconfs=clust_size*njacks;
  
  cout<<"NTimes: "<<ntimes<<endl;
  cout<<"NFlavs: "<<nflavs<<endl;
  cout<<"Nops: "<<nops<<endl;
  cout<<"nconf_types: "<<nconf_types<<endl;
  cout<<"NCopies: "<<ncopies<<endl;
  cout<<"NConfs_max: "<<nconfs_max<<endl;
  cout<<"Clust_size: "<<clust_size<<endl;
  cout<<"Nconfs: "<<nconfs<<endl;
  
  const size_t sigma_xy_tag=6,i_sigma_xy=get_op_pos(op_list,sigma_xy_tag);
  const size_t sigma_zt_tag=9,i_sigma_zt=get_op_pos(op_list,sigma_zt_tag);
  
  const index_t ind({{"conf",nconfs},{"copy",ncopies},{"tmeas",ntimes},{"flav",nflavs},{"op",nops},{"conf_type",nconf_types}});
  
  cout<<"Operator Sigma_xy pos: "<<i_sigma_xy<<endl;
  cout<<"Operator Sigma_zt pos: "<<i_sigma_zt<<endl;
  
  enum{UP,DW,ST};
  
  index_t ind_res({{"time",ntimes},{"conf_type",nconf_types},{"flav",nflavs}});
  djvec_t C_f(ind_res.max());
  
  auto time=tmeas_list.begin();
  for(size_t itime=0;itime<ntimes;itime++)
    for(size_t conf_type=0;conf_type<nconf_types;conf_type++)
      {
	string conf_type_tag[2]={"gauge","ferm"};
	
	cout<<endl;
	cout<<"======================================================="<<endl;
	cout<<"time: "<<(*time)<<endl;
	cout<<"conf_type: "<<conf_type_tag[conf_type]<<endl;
	
	djack_t plaq=get_obs(obs_list,ind,itime,/*iflav*/0,i_sigma_xy,conf_type,PLAQ);
	cout<<"Plaq: "<<smart_print(plaq)<<endl;
	
	djack_t Q=get_obs(obs_list,ind,itime,/*iflav*/0,i_sigma_xy,conf_type,TOT_CHARGE);
	cout<<"Q: "<<smart_print(Q)<<endl;
	
	for(size_t iflav=0;iflav<nflavs;iflav++)
	  {
	    cout<<endl;
	    cout<<"-------------------------------------------------------"<<endl;
	    cout<<"iflav: "<<iflav<<endl;
	    
	    djack_t sigma_xy=get_obs(obs_list,ind,itime,iflav,i_sigma_xy,conf_type,TENS_IM);
	    djack_t Qsigma_zt=get_obs(obs_list,ind,itime,iflav,i_sigma_zt,conf_type,SPI_IM);
	    djack_t rat=Qsigma_zt/sigma_xy*1000;
	    cout<<"rat: "<<smart_print(rat)<<endl;
	    
	    const size_t glb_vol=L*L*L*T;
	    const djack_t Q2=get_obs(obs_list,ind,itime,iflav,i_sigma_xy,conf_type,TOT_CHARGE2)/glb_vol;
	    cout<<"Q2: "<<smart_print(Q2)<<endl;
	    
	    size_t ires=ind_res({itime,conf_type,iflav});
	    C_f[ires]=Qsigma_zt/(sigma_xy*sqrt(Q2));
	    
	    cout<<"C_f: "<<smart_print(C_f[ires])<<endl;
	  }
	
	time++;
      }
  
  //write the plot
  grace_file_t out_C_f("C_f.xmg");
    for(size_t conf_type=0;conf_type<nconf_types;conf_type++)
      {
	out_C_f.new_data_set();
	for(size_t itime=0;itime<ntimes;itime++)
	  out_C_f.write_ave_err(itime,C_f[ind_res({itime,conf_type,0})].ave_err());
      }
    
  return 0;
}
