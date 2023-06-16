#include <tranalisi.hpp>

const int T=128,L=64,TH=T/2,tw=22;
const range_t file_range={500,8,2100};

enum{P5P5,S0P5,V0P5,V1P5,V2P5,V3P5,A0P5,A1P5,A2P5,A3P5};
enum{VNU1,VNU2,ANU0,ANU1,ANU2,ANU3};
enum{VNUIJ,ANU30,ANUII,ANU03,ANU33};
enum{CNS,LOC};
const size_t nContrs=ANU33+1;

std::vector<std::string> thTag{"0.00000","-1.637960","1.637960","-4.106540","4.106540","-5.767930","5.767930"};
std::vector<double> thList;
const size_t nth=thTag.size(),nm=2;
const size_t billOfMu[]{V0P5,V1P5,V2P5,V3P5};

int main()
{
  for(const auto& thT : thTag)
    thList.push_back(strtod(thT.c_str(),nullptr));
  
  set_njacks(50);
  
  index_t id3pts({{"m",2},{"th",nth},{"cnsLoc",2},{"contr",nContrs}});
  
  std::vector<djvec_t> c3pts(id3pts.max(),djvec_t(T));

  const string dataPath{"c3pts.dat"};
  raw_file_t dataFile;
  
  if(file_exists(dataPath))
    {
      dataFile.open(dataPath,"r");
      for(auto& c : c3pts)
	dataFile.bin_read(c);
    }
  else
    {
      for(size_t ith=0;ith<nth;ith++)
	{
	  auto get3ptGetter=[ith](const string& cnsLoc,
				  const int& nmu)
	  {
	    return
	      [data=read_conf_set_t("out/%04d/mes_contr_PHOTON_EMISSION_TH_"+thTag[ith]+"_M0_0.0184_M1_0.237"+cnsLoc,file_range,2,{0,1},T,VERBOSE),
	       id=index_t({{"m1",2},{"avnu",6},{"mu",nmu},{"bil",10},{"reim",2},{"t",T}})]
	      (const size_t& im,
	       const size_t& avnu,
	       const size_t& mu,
	       const size_t& bil,
	       const size_t& reim)
	      {
		djvec_t res(T);
		for(size_t t=0;t<T;t++)
		  res[t]=data[id({im,avnu,mu,bil,reim,t})];
		
		return res;
	      };
	  };
	  
	  auto getLoc=get3ptGetter("_LOC",1);
	  auto getCns=get3ptGetter("",4);
	  
	  for(size_t im=0;im<nm;im++)
	    {
	      for(const auto& [contr,avnu,mu,reim,f] :
		    {std::make_tuple
		     (VNUIJ,VNU1,2,IM,0.5),
		     {VNUIJ,VNU2,1,IM,-0.5},
		     {ANU30,ANU3,0,IM,1.0},
		     {ANUII,ANU1,1,RE,0.5},
		     {ANUII,ANU2,2,RE,0.5},
		     {ANU03,ANU0,3,IM,1.0},
		     {ANU33,ANU3,3,RE,1.0}})
		{
		  c3pts[id3pts({im,ith,CNS,contr})]+=f*getCns(im,avnu,mu,S0P5,reim);
		  c3pts[id3pts({im,ith,LOC,contr})]+=f*getLoc(im,avnu,0,billOfMu[mu],reim);
		}
	    }
	}
      
      dataFile.open(dataPath,"w");
      for(auto& c : c3pts)
	dataFile.bin_write(c);
    }
  
  djvec_t fP(4);
  for(size_t cnsLoc : {CNS,LOC})
    for(size_t im=0;im<2;im++)
      {
	djack_t& _fP=fP[cnsLoc+2*im];
	_fP=0.0;
	c3pts[id3pts({im,0/*th*/,cnsLoc,ANUII})].ave_err().write("plots/fP_corr_im_"+to_string(im)+"_consloc_"+to_string(cnsLoc)+".xmg");
	for(size_t t=0;t<T;t++)
	  _fP+=c3pts[id3pts({im,0/*th*/,cnsLoc,ANUII})][t];
	
	cout<<"fP[im="<<im<<",cnsLoc="<<cnsLoc<<"]: "<<_fP.ave_err()<<endl;
      }
  
  for(size_t im=0;im<2;im++)
    {
      const djack_t Zv=fP[CNS+2*im]/fP[LOC+2*im];
      cout<<"Zv[im="<<im<<"]: "<<Zv.ave_err()<<endl;
    }
  
  for(size_t i=0;i<id3pts.max();i++)
    {
      const vector<size_t> comps=id3pts(i);
      const size_t im=comps[0],cnsLoc=comps[2];
      c3pts[i]/=fP[cnsLoc+2*im];
      c3pts[i].ave_err().write("plots/"+id3pts.descr(i)+".xmg");
      
      djvec_t c=c3pts[i].subset(tw,tw+TH);
      const djack_t E=constant_fit(effective_mass(c),20,30,"plots/eff_mass_"+id3pts.descr(i)+".xmg");
      
      for(size_t it=0;it<=TH;it++)
	c[it]/=exp(-E.ave()*it);
      
      c.ave_err().write("plots/prec_"+id3pts.descr(i)+".xmg");
    }
  
  const double a=0.4;
  const double eps=0.4*a;
  const double M=0.81;
  const size_t im=1;
  for(size_t ith=0;ith<nth;ith++)
    {
      const double th=thList[ith];
      const double k=M_PI*th/L;
      const double kHat=2*sin(k/2);
      const double eg_v0=2*asinh(fabs(kHat)/2);
      const double xg=2*eg_v0/M;
      
      cout<<"ith: "<<ith<<", th: "<<th<<", xg: "<<xg<<endl;
      const size_t nxk=39;
      for(size_t ixk=0;ixk<=nxk;ixk++)
	{
	  const double xk=ixk/(double)nxk;
	  const double virt=M*xk;
	  const double eg=sqrt(sqr(eg_v0)+sqr(virt));
	  
	  grace_file_t re("plots/re"+id3pts.descr(id3pts({im,ith,LOC,VNUIJ}))+"_ixk"+std::to_string(ixk)+".xmg");
	  re.write_line([eg,eps](const double e){return sinh(e-eg)/(sqr(sinh(e-eg)+sqr(eps)));},0.0,1,grace::RED);
	}
    }
  
  //1/(sinh(x-eg)-i*eps);
  
      // effective_mass(am0nu11,T/2,0).ave_err().write("plots/m0_a11.xmg");
      // const djvec_t am1nu11=get3pt(1,ANU1,1,S0P5,RE);
      // effective_mass(am1nu11,T/2,0).ave_err().write("plots/m1_a11.xmg");
      
      // auto get2pt=
      //   [data=read_conf_set_t("out/%04d/mes_contr_MES_2PT_NO_SM",file_range,2,{0},T,VERBOSE),
      //    id=index_t({{"m01",4},{"bil",10},{"t",T}})]
      //   (const size_t& im01,
      //    const size_t& bil)
  //   {
  //     djvec_t res(T);
  //     for(size_t t=0;t<T;t++)
  // 	res[t]=data[id({im01,bil,t})];
      
  //     return res;
  //   };
  
  // const djvec_t Pcs=get2pt(0,P5P5).symmetrized();
  // //const djvec_t Pss=get2pt(1,P5P5).symmetrized();
  // const djvec_t Pss=get2pt(2,P5P5).symmetrized();
  // const djvec_t Pcc=get2pt(3,P5P5).symmetrized();
  // effective_mass(Pcs).ave_err().write("plots/cPcs.xmg");
  // effective_mass(Pss).ave_err().write("plots/cPss.xmg");
  // effective_mass(Pcc).ave_err().write("plots/cPcc.xmg");
  
  return 0;
}
