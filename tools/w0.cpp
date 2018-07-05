#include <tranalisi.hpp>

#include <set>

djvec_t load(const string &path,const size_t ncols)
{
  raw_file_t fin(path,"r");
  
  set<size_t> tlist;
  vector<double> E;
  char line[1024];
  while(fin.get_line(line)!=NULL)
    {
      size_t iconf,imeas,acc;
      double p,r;
      double e;
      switch(ncols)
	{
	case 6:
	  if(sscanf(line,"%zu %zu %zu %lg %lg %lg",&iconf,&acc,&imeas,&p,&r,&e)!=6) CRASH("parsing %s",line);
	  break;
	case 4:
	  if(sscanf(line,"%zu %zu %zu %lg",&iconf,&acc,&imeas,&e)!=4) CRASH("parsing %s",line);
	  break;
	default:
	  CRASH("Unknwon");
	}
      tlist.insert(imeas);
      //const double c1=-1.0/12,c0=1-8*c1;
      //E.push_back(c0*6.0*(1.0-p)+c1*12.0*(1.0-r));
      E.push_back(e);
    };
  fin.close();
  
  const size_t nt=tlist.size();
  const size_t nlines=E.size();
  const size_t nconfs_poss=nlines/nt;
  const size_t clust_size=nconfs_poss/njacks;
  const size_t nconfs=clust_size*njacks;
  
  if(nconfs_poss*nt!=nlines) CRASH("file is %zu line length, not a multiple of nt %zu",nlines,nt);
  
  cout<<"NLines: "<<nlines<<endl;
  cout<<"Nt: "<<nt<<endl;
  cout<<"NConfsPoss: "<<nconfs_poss<<endl;
  cout<<"NJacks: "<<njacks<<endl;
  cout<<"ClustSize: "<<clust_size<<endl;
  cout<<"NConfs: "<<nconfs<<endl;
  if(nconfs==0) CRASH("Needs at least %zu confs",njacks);
  
  //put int jacks
  djvec_t res(nt);
  res=0.0;
  
  index_t conf_t_id({{"Conf",nconfs_poss},{"T",nt}});
  for(size_t iconf=0;iconf<nconfs;iconf++)
    {
      const size_t iclust=iconf/clust_size;
      for(size_t it=0;it<nt;it++)
	res[it][iclust]+=E[conf_t_id({iconf,it})];
    }
  res.clusterize(clust_size);
  
  return res;
}

djack_t find_w0(const vector<double> &t,const djvec_t &E,const size_t ord=5,const double w0_ref=0.3)
{
  const double dt=t[1]-t[0];
  const size_t nt=E.size();
  
  djvec_t t2E(nt);
  for(size_t it=0;it<nt;it++) t2E[it]=E[it]*sqr(t[it]);
  t2E.ave_err().write(t,"plots/t2E.xmg");
  
  //compute derivative of t2^p2 and multiply by t
  djvec_t r(nt);
  for(size_t it=0;it<nt;it++)
    {
      bool bw=(it>0),fw=(it+1<nt);
      djack_t d=0.0;
      if(bw) d+=(t2E[it]-t2E[it-1])/dt;
      if(fw) d+=(t2E[it+1]-t2E[it])/dt;
      d/=(bw+fw);
      
      r[it]=d*t[it];
    }
  r.ave_err().write(t,"plots/tdt2E.xmg");
  
  double guess_t=0.0;
  //finds where average pass ref
  for(size_t it=0;it<t.size();it++)
    if(r[it].ave()<w0_ref) guess_t=t[it]+dt*0.5;
  
  if(ord%2!=1) CRASH("only odd ords implemented");
  const double dfit=(ord+1)*dt*0.5;
  const double guess_t_min=guess_t-dfit;
  const double guess_t_max=guess_t+dfit;
  djvec_t coeffs=poly_fit(t,r,ord,guess_t_min,guess_t_max);
  
  djack_t w02_fr_a2;
  for(size_t ijack=0;ijack<=njacks;ijack++)
    w02_fr_a2[ijack]=Brent_solve(
				 [w0_ref,&coeffs,ijack](double x)
				 {
				   double out=-w0_ref,e=1;
				   for(size_t deg=0;deg<coeffs.size();deg++)
				     {
				       out+=coeffs[deg][ijack]*e;
				       e*=x;
				     }
				   return out;
				 },
				 guess_t_min,guess_t_max);
  
  cout<<"Intercept: "<<smart_print(w02_fr_a2.ave_err())<<endl;
  
  //print the fit
  grace_file_t fout("plots/fit.xmg");
  fout.write_polygon(bind(poly_eval<djvec_t>,coeffs,_1),guess_t_min,guess_t_max);
  fout.write_line([w0_ref](double x){return w0_ref;},guess_t_min,guess_t_max);
  fout.write_vec_ave_err(t,r.ave_err());
  fout.set_settype(grace::XYDX);
  fout<<w02_fr_a2.ave()<<" "<<w0_ref<<" "<<w02_fr_a2.err()<<"\n";
  
  return w02_fr_a2;
}

int main(int narg,char **arg)
{
  if(narg<5) CRASH("Use %s path dt njacks ncols",arg[0]);
  string path=arg[1];
  double dt=strtod(arg[2],NULL);
  set_njacks(atoi(arg[3]));
  size_t ncols=atoi(arg[4]);
  
  const djvec_t p=load(path,ncols);
  const size_t nt=p.size();
  
  //compute times, t^2*p
  vector<double> t(nt);
  djvec_t t2p(nt);
  for(size_t it=0;it<nt;it++)
    {
      t[it]=dt*it;
      t2p[it]=sqr(t[it])*p[it];
    }
  p.ave_err().write(t,"plots/plaq.xmg");
  
  const djack_t w02_fr_a2=find_w0(t,p);
  const double w0=0.1755;
  djack_t a=w0/sqrt(w02_fr_a2);
  cout<<"a: "<<smart_print(a)<<" fm"<<endl;
  
  return 0;
}
