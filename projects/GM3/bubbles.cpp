#include "index.hpp"
#include "jack.hpp"
#include "meas_vec.hpp"
#include "raw_file.hpp"

const size_t T=224;
const size_t nDirs=3;
const size_t nHits=10;
const size_t nCopies=8;
const size_t nConfs=10;
const char confs[nConfs][11]={"0008_r1","1100_r0","1108_r0","1132_r0","1148_r0","1180_r0","1188_r0","1196_r0","1204_r0","1212_r0"};

djvec_t load(const std::vector<std::pair<const char*,const char*>>& list)
{
  std::map<const char*,double> mass{{"l",0.00044},{"l2",0.00088},{"l3",0.00132},{"l4",0.00220},{"s",0.011787}};
  
  index_t dataId({{"conf",nConfs},{"t",T}});
  vector<double> confData(dataId.max());
  
  for(size_t iConf=0;iConf<nConfs;iConf++)
    {
      index_t loadId({{"dir",nDirs},{"cj",nCopies*nHits},{"t",T}});
      vector<double> d(loadId.max());
	for(size_t r=0;r<2;r++)
	  for(const auto& [q1,q2] : list)
	    {
	      raw_file_t in(combine("/home/francesco/QCD/LAVORI/GM4/bubbles/E112/data/%s/mes_contr_OET_%s%s_R%zu",confs[iConf],q1,q2,r),"r");
	      line_t line;
	      
	      auto getLine=
		[&in,
		 &line]()
		{
		  if(not in.get_line(line))
		    CRASH("error getting line on file %s",in.get_path().c_str());
		};
	      
	      for(size_t icopy=0;icopy<nCopies;icopy++)
		for(size_t ihit=0;ihit<nHits;ihit++)
		  {
		    getLine();
		    size_t t0,x0,y0,z0;
		    char firstProp[54];
		    char secondProp[54];
		    if(const int n=sscanf(line," # Contraction of %s ^ \\dag and %s origin located at txyz = (%zd,%zd,%zd,%zd)",firstProp,secondProp,&t0,&x0,&y0,&z0);n!=6)
		      CRASH("%d %s",n,line);
		    printf("%s %s %zu %zu %zu %zu\n",firstProp,secondProp,t0,x0,y0,z0);
		    
		    bool ab[2];
		    const char* q[2]={q1,q2};
		    for(int ip=0;ip<2;ip++)
		      ab[ip]=
			combine("S_m%s_R%zu_copy%zu",q[ip],r,icopy)==firstProp and
			combine("S_m%s_R%zu_copy%zu,",q[1-ip],r,icopy)==secondProp;
		    
		    if(not (ab[0] xor ab[1]))
		      CRASH("%d %d     %s %s",ab[0],ab[1],firstProp,secondProp);
		    // printf("%s %s %d %d\n",q1,q2,ab[0],ab[1]);
		    
		    //getLine();
		    
		    for(size_t dir=0;dir<3;dir++)
		      {
			getLine();
			size_t mu;
			if(const int n=sscanf(line," # A%zdP5",&mu);n!=1)
			  CRASH("%d %s",n,line);
			if(mu!=dir+1)
			  CRASH("%s %zu",line,mu);
			
			for(size_t t=0;t<T;t++)
			  {
			    getLine();
			    double im;
			    if(const int n=sscanf(line,"%lg",&im);n!=1)
			      CRASH("%d %s",n,line);
			    d[loadId({dir,ihit+nHits*icopy,(t+t0)%T})]+=im/2*(mass[q1]+mass[q2])*(ab[1]?-1:1);
			  }
			// getLine();
		      }
		  }
	    }
	
	// vector<double> sum(T,0.0);
	// for(size_t dir=0;dir<nDirs;dir++)
	//   for(size_t i=0;i<nCopies*nHits;i++)
	//     for(size_t j=0;j<nCopies*nHits;j++)
	// 	if(i!=j)
	// 	  {
	// 	    for(size_t t1=0;t1<T;t1++)
	// 	      for(size_t t2=0;t2<T;t2++)
	// 		{
	// 		  size_t dt=(t2-t1+T)%T;
	// 		  sum[dt]+=d[id({dir,i,t1})]*d[id({dir,j,t2})];
	// 		}
	// 	    printf("%zu %zu %zu\n",dir,i,j);
	// 	  }
	
	// for(size_t dt=0;dt<T;dt++)
	//   printf("%zu %lg\n",dt,sum[dt]);
	
	vector<double> sum(T,0.0);
	for(size_t dir=0;dir<nDirs;dir++)
	  {
	    vector<double> A(T,0.0);
	    vector<double> AA(T,0.0);
	    for(size_t i=0;i<nCopies*nHits;i++)
	      for(size_t t=0;t<T;t++)
		{
		  A[t]+=d[loadId({dir,i,t})];
		  for(size_t dt=0;dt<T;dt++)
		    AA[dt]+=d[loadId({dir,i,t})]*d[loadId({dir,i,(t+dt)%T})];
		}
	    
	    for(size_t dt=0;dt<T;dt++)
	      {
		for(size_t t=0;t<T;t++)
		  sum[dt]+=A[t]*A[(t+dt)%T];
		sum[dt]-=AA[dt];
	      }
	  }
	
	for(size_t dt=0;dt<T;dt++)
	  confData[dataId({iConf,dt})]=sum[dt];
    }
  
  djvec_t data(T);
  jackknivesFill(nConfs,
		 [&data,
		  &dataId,
		  &confData](const size_t& iConf,
			     const size_t& iClust,
			     const double& w)
		 {
		   for(size_t dt=0;dt<T;dt++)
		     data[dt][iClust]+=confData[dataId({iConf,dt})];
		 });
  data.clusterize((double)nConfs/njacks);
  
  return data;
}

int main()
{
  set_njacks(nConfs);
  
  const djvec_t l234s=load({{"l","l2"},{"l2","l3"},{"l3","l4"},{"l4","s"}});
  const djvec_t ls=load({{"l","s"}});
  l234s.ave_err().write("/tmp/l234s.xmg");
  ls.ave_err().write("/tmp/ls.xmg");
  // const djvec_t l3l4=load("l3","l4");
  // const djvec_t l4s=load("l4","s");
  // const djvec_t ls=load("l","s");
  
  // (ll2+l2l3+l3l4+l4s).ave_err().write("/tmp/ls_ps.xmg");
  
  return 0;
}
