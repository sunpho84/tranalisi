#include "ave_err.hpp"
#include "common.hpp"
// #include "gevp.hpp"

bool insOnDs;

size_t L;

size_t tPi;
size_t tpExpMin;
size_t tpExpMax;
size_t locNjacks;
vector<size_t> twList;
double mFinalMes;
vector<size_t> LS_fit_range_min(3);
vector<size_t> LS_fit_range_max(3);

const int nls=1;
double ml;
double ms;
size_t nmc;
vector<double> mcList;

const string plotsPath="plots";

int main()
{
  {
    input_file_t input("input.txt");
    input.read(L,"L");
    input.read(T,"T");
    input.read(tPi,"tPi");
    input.read(tpExpMin,"tpExpMin");
    input.read(tpExpMax,"tpExpMax");
    input.expect("LSFitRange");
    for(size_t i=0;i<3;i++)
      input.read(LS_fit_range_min[i]);
    for(size_t i=0;i<3;i++)
      input.read(LS_fit_range_max[i]);
    
    input.read(locNjacks,"NJacks");
    input.read(insOnDs,"insOnDs");
    input.read(ml,"ml");
    input.read(ms,"ms");
    input.read(mFinalMes,"mFinalMes");
    nmc=input.read<size_t>("nmc");
    mcList.resize(nmc);
    for(size_t ic=0;ic<nmc;ic++)
      mcList[ic]=input.read<double>();
    
    size_t locNtw;
    input.read(locNtw,"ntw");
    twList.resize(locNtw);
    for(size_t itw=0;itw<locNtw;itw++)
      input.read(twList[itw]);
  }
  
  setConfs();
  set_njacks(locNjacks);
  cout<<"NConfs: "<<nConfs<<endl;
  
  map<string,array<djvec_t,2>> data;
  for(const char* entry : {"2p","3p","2pv","3pUnsme"
    })
    {
      const auto rawData=getRawData(entry);
      
      for(const auto& [tag,v] : rawData)
	{
	  grace_file_t f(plotsPath+"/full_"+tag+".xmg");
	  for(size_t iConf=0;iConf<nConfs;iConf++)
	    {
	      for(size_t t=0;t<T;t++)
		f.write_xy(t,v[iConf][t]);
	      f.set_legend(confs[iConf]);
	      f.new_data_set();
	    }
	  f.close();
	}
      
      for(const auto& it : rawData)
	for(size_t ri=0;ri<2;ri++)
	  {
	    const string& tag=it.first;
	    const vector<vector<double>>& v=it.second;
	    djvec_t& d=data[tag][ri];
	    d.resize(T);
	    jackknivesFill(nConfs,[&v,
				   &ri,
				   &d](const size_t& iConf,
				       const size_t& iClust,
				       const double& weight)
	    {
	      for(size_t t=0;t<T;t++)
		d[t][iClust]+=v[iConf][t+T*ri]*weight;
	    });
	  d.clusterize((double)nConfs/njacks);
	}
    }
  
  const map<string,ave_err_t> liter=
    {{"fnalpD0ToPi",{3.119,0.057}},
     {"fnal0D0ToPi",{1.2783,0.0061}},
     {"etmpD0ToPi",{2.130,0.096}},
     {"etm0D0ToPi",{1.134,0.049}},
     {"fnalpD0ToK",{1.451,0.017}},
     {"fnal0D0ToK",{1.0240,0.0021}},
     {"etmpD0ToK",{1.336,0.054}},
     {"etm0D0ToK",{0.979,0.019}},
     {"fnalpD0sToK",{1.576,0.013}},
     {"fnal0D0sToK",{0.9843,0.0018}}};
  
  auto get=
    [&data](const string& name,
	    const bool reim) -> auto&
    {
      auto f=data.find(name);
      if(f==data.end())
	CRASH("unable to find %s",name.c_str());
      
      return f->second[reim];
    };
  
  auto getP5P5=
    [&get](const string& name,
	    const bool reim) -> auto&
    {
      return get(name+"__P5P5",reim);
    };
  
  auto getVKVK=
    [&get](const string& name,
	    const bool reim)
    {
      return (get(name+"__V1V1",reim)+
	      get(name+"__V2V2",reim)+
	      get(name+"__V3V3",reim))/3.0;
    };
  
  djvec_t ZPLList(nls),MPList(nls);
  for(size_t ils=0;ils<nls;ils++)
    {
      const string ls=(ils==0)?"l":"s";
      const djvec_t P=getP5P5("S_m"+ls+"__S_ml",0).symmetrized();
      djack_t ZPLL;
      const string Pname=ils?"K":"Pi";
      two_pts_fit(ZPLL,MPList[ils],P,T/2,tPi,tPi+20,plotsPath+"/"+Pname+".xmg");
      ZPLList[ils]=sqrt(ZPLL);
      cout<<"MP"+ls+":     "<<smart_print(MPList[ils])<<endl;
      
      // const djvec_t Pdt=getP5P5("S_m"+ls+"_dth__S_ml",1).symmetrized();
      // Pdt.ave_err().write(plotsPath+"/"+Pname+"dth.xmg");
    }
  
  const djack_t fP=ZPLList[0]*2*ml/sqr(MPList[0]);
  
  const djack_t a=MPList[0]/mFinalMes;
  cout<<"a: "<<a<<endl;
  cout<<"fpi: "<<(fP/a).ave_err()<<" GeV"<<endl;
  
  const djack_t xi=MPList[0]/fP;
  cout<<"xi: "<<xi.ave_err()<<endl;
  
  // if(not file_exists("input"))
  //   CRASH("input file not present");
  // input_file_t input("input");
  // const size_t tw=input.read<size_t>("tw");
  const size_t ntw=twList.size();
  const size_t nsme=1;
  const size_t isme=0;
  
  std::vector<pair<size_t,size_t>> lsCombo{{0,0}};
  if(nls>1)
    for(const auto& [l,s] : {std::make_pair(0,1),{1,0}})
      lsCombo.emplace_back(l,s);
  const size_t nlsCombo=lsCombo.size();
  
  const index_t iDs({{"ls",nls},{"c",nmc},{"sme",nsme}});
  const index_t iEstim({{"tw",ntw},{"lsCombo",nlsCombo},{"c",nmc},{"sme",nsme}});
  
  auto getObs=
    [&iDs]()
    {
      return
	[&iDs,
	 data=djvec_t(iDs.max())](const size_t& ils,
				  const size_t& ic,
				  const size_t& isme) mutable -> djack_t&
	{
	  return data[iDs({ils,ic,isme})];
	};
    };
  
  auto getEstim=
    [&iEstim]()
    {
      return
	[&iEstim,
	 data=std::vector<djvec_t>(iEstim.max(),djvec_t())](const size_t& itw,
				     const size_t& ilsCombo,
				     const size_t& ic,
				     const size_t& isme) mutable -> djvec_t&
	{
	  return data[iEstim({itw,ilsCombo,ic,isme})];
	};
    };
  
  auto MD=getObs();
  auto ZDL=getObs();
  auto ZDS=getObs();
  auto delta=getEstim();
  auto cfpList=getEstim();
  auto cf0List=getEstim();
  
  for(size_t ils=0;ils<nls;ils++)
    for(size_t ic=0;ic<nmc;ic++)
      {
	const string Dname="D"+to_string(ic)+((ils==1)?"s":"");
	
	const string ls=((ils==0)?"l":"s");
	
	djack_t M0;
	const djvec_t D=getP5P5("S_m"+ls+"__S_mc"+to_string(ic),0).symmetrized();
	const djvec_t Dstar=getVKVK("S_m"+ls+"__S_mc"+to_string(ic),0).symmetrized();
	effective_mass(Dstar).ave_err().write(plotsPath+"/"+Dname+"star.xmg");
	
	djack_t ZDLL0;
	two_pts_fit(ZDLL0,M0,D,T/2,27,30,plotsPath+"/"+Dname+".xmg");
	
	djvec_t excD=D;
	for(size_t t=0;t<excD.size();t++)
	  excD[t]-=two_pts_corr_fun(ZDLL0,M0,T/2.0,t,1);
	
	djack_t ZDLL1,M1;
	two_pts_fit(ZDLL1,M1,excD,T/2,14,27,plotsPath+"/exc_"+Dname+".xmg");
	
	jack_fit_t fitter;
	djvec_t pars(4);
	fitter.add_fit_par_limits(pars[0],"Z",ZDLL0.ave(),ZDLL0.err(),ZDLL0.ave()*0.9,ZDLL0.ave()*1.1);
	fitter.add_fit_par_limits(pars[1],"M",M0.ave(),M0.err(),M0.ave()*0.9,M0.ave()*1.1);
	fitter.add_fit_par_limits(pars[2],"Z1",ZDLL1.ave(),ZDLL1.err(),fabs(ZDLL1.ave())*0.9,fabs(ZDLL1.ave())*1.1);
	fitter.add_fit_par_limits(pars[3],"M1",M1.ave(),M1.err(),M1.ave()*0.9,M1.ave()*1.6);
	for(size_t it=14;it<=30;it++)
	  fitter.add_point(D[it],[t=(double)it](const vector<double>& pars,
						  const size_t iel)
	  {
	    return two_pts_corr_fun(pars[0],pars[1],T/2.0,t,1)+two_pts_corr_fun(pars[2],pars[3],T/2.0,t,1);
	  });
	fitter.fit();
	for(const auto& p : pars)
	  cout<<""<<p.ave_err()<<endl;
	
	MD(ils,ic,0)=pars[1];
	
	for(size_t is=0;is<nsme;is++)
	  {
	    const djvec_t D_sme=getP5P5("S_m"+ls+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	    const djvec_t D_sme_sme=getP5P5("S_m"+ls+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	    
	    djack_t M;
	    djack_t ZDLS,ZDSS;
	    two_pts_fit(ZDLS,M,D_sme,T/2,22,28,plotsPath+"/"+Dname+"_sme"+to_string(is)+".xmg");
	    two_pts_fit(ZDSS,M,D_sme_sme,T/2,22,28,plotsPath+"/"+Dname+"_sme"+to_string(is)+"_sme"+to_string(is)+".xmg");
	    
	    djack_t& zl=ZDL(ils,ic,is);
	    djack_t& zs=ZDS(ils,ic,is);
	    
	    /// LS FIT
	    {
	      //parameters to fit
	      minimizer_pars_t pars;
	      pars.add("M",M.ave(),M.err());
	      pars.add("Z_L",sqrt(ZDLL0).ave(),sqrt(ZDLL0).err());
	      pars.add("Z_S",sqrt(ZDSS).ave(),sqrt(ZDSS).err());
	      
	      //! fit for real
	      size_t iel=0;
	      auto x=vector_up_to<double>(D.size());
	      multi_ch2_t<djvec_t> two_pts_LS_fit_obj({x,x,x},LS_fit_range_min,LS_fit_range_max,{D,D_sme,D_sme_sme},
						      {two_pts_corr_fun_t(T/2,1),two_pts_corr_fun_t(T/2,1),two_pts_corr_fun_t(T/2,1)},
						      [](const vector<double> &p,size_t icontr)
						      {
							switch(icontr)
							  {
							  case 0:return vector<double>({p[1]*p[1],p[0]});break;
							  case 1:return vector<double>({p[1]*p[2],p[0]});break;
							  case 2:return vector<double>({p[2]*p[2],p[0]});break;
							  default: CRASH("Unknown contr %zu",icontr);return p;
							  }
						      },iel);
	      
	      minimizer_t minimizer(two_pts_LS_fit_obj,pars);
	      
	      for(iel=0;iel<njacks;iel++)
		{
		  //minimize and print the result
		  vector<double> par_min=minimizer.minimize();
		  M[iel]=par_min[0];
		  zl[iel]=par_min[1];
		  zs[iel]=par_min[2];
		}
	      
	      grace_file_t out(plotsPath+"/"+Dname+"_SLfit_"+to_string(is)+"_sme"+to_string(is)+".xmg");
	      out.write_polygon([&M](double x){return M;},18,23);
	      out.write_vec_ave_err(x,effective_mass(D,T/2,1).ave_err());
	      out.write_vec_ave_err(x,effective_mass(D_sme,T/2,1).ave_err());
	      out.write_vec_ave_err(x,effective_mass(D_sme_sme,T/2,1).ave_err());
	    }
	    
	    cout<<"M"<<Dname<<": "<<M.ave_err()<<" "<<(M/a).ave_err()<<" (lattice self-estimated)"<<endl;
	    cout<<"ZDL: "<<zs.ave_err()<<endl;
	    cout<<"ZDS: "<<zl.ave_err()<<endl;
	    cout<<"ZDLS: "<<ZDLS.ave_err()<<" ZDL*ZDS: "<<(zs*zl).ave_err()<<endl;
	    
	    const djvec_t Dstar_sme=getVKVK("S_m"+ls+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	    const djvec_t Dstar_sme_sme=getVKVK("S_m"+ls+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
	    effective_mass(Dstar_sme).ave_err().write(plotsPath+"/"+Dname+"star_sme"+to_string(is)+".xmg");
	    const djack_t MDStar=constant_fit(effective_mass(Dstar_sme_sme),22,28,plotsPath+"/"+Dname+"star_sme"+to_string(is)+"_sme"+to_string(is)+".xmg");
	    cout<<"MDStar: "<<MDStar.ave_err()<<endl;
	    
	    // const size_t t0=3;
	    // auto [eig,recastEigvec,origEigvec]=
	    //   gevp({D,D_sme,D_sme,D_sme_sme},t0);
	    
	    // effective_mass(eig[0]).ave_err().write(plotsPath+"/"+Dname+"_eig1.xmg");
	    // effective_mass(eig[1]).ave_err().write(plotsPath+"/"+Dname+"_eig2.xmg");
	    
	    const djvec_t& k=delta(0,ils,ic,isme)=effective_mass(Dstar_sme_sme)-M+MPList[0];
	    k.ave_err().write(plotsPath+"/"+Dname+"delta.xmg");
	  }
	
      }
  
  if(nsme>1)
    CRASH("Too many sme");
  
  const size_t dT=60*0.40/a.ave();
  
  for(size_t ilsCombo=0;ilsCombo<lsCombo.size();ilsCombo++)
    {
      const size_t& ilsSeq=lsCombo[ilsCombo].first;
      const size_t& ilsSpe=lsCombo[ilsCombo].second;
      
      const string lsSeq=((ilsSeq==0)?"l":"s");
      const string lsSpe=((ilsSpe==0)?"l":"s");
      
      const size_t i=ilsSeq+ilsSpe;
      const djack_t MP=MPList[i];
      const djack_t ZP=ZPLList[ilsSeq+ilsSpe];
      const string Pname=(ilsSpe==1 or ilsSeq==1)?"K":"Pi";
      
      const double mseq=(ilsSeq==0)?ml:ms;
      
      const size_t tmin3=27,tmax3=35;
      djvec_t f0(nmc),fp(nmc),ft(nmc);
      for(size_t ic=0;ic<nmc;ic++)
	{
	  const string Dname="D"+to_string(ic)+((ilsSpe==1)?"s":"");
	  const double mc=mcList[ic];
	  
	  djack_t& M=MD(ilsSpe,ic,0);
	  const djack_t Q2Max=sqr(M-MP)/sqr(a);
	  cout<<Dname<<" -> "<<Pname<<" Q2 max: "<<Q2Max.ave_err()<<endl;
	  
	  for(size_t is=0;is<nsme;is++)
	    {
	      // const djvec_t D_sme=getP5P5("S_m"+lsSpe+"_sme"+to_string(is)+"__S_mc"+to_string(ic),0).symmetrized();
	      
	      // const djvec_t N=D_sme/two_pts_corr_fun(ZDL(ilsSpe,ic,is)*ZDS(ilsSpe,ic,is),M,T/2.,0,1);
	      for(size_t itw=0;itw<ntw;itw++)
		{
		  const size_t tw=twList[itw];
		  auto g=
		    [&](const string& name,
			const bool ri,
			auto& ZD)
		    {
		      const djvec_t& R=getP5P5(name,ri);
		      djvec_t s(T);
		      for(size_t t=0;t<T;t++)
			s[t]=R[(t+tw)%T];
		      s.ave_err().write(plotsPath+"/raw_"+name+".xmg");
		      
		      const djvec_t e=(-symmetric_derivative(log(s)));
		      grace_file_t eff_plot(plotsPath+"/eff_"+name+".xmg");
		      eff_plot.set_no_line();
		      eff_plot.write_vec_ave_err(e.ave_err());
		      eff_plot.write_vec_ave_err((-e).symmetric().ave_err());
		      eff_plot.set_all_colors(grace::BLUE);
		      
		      djvec_t A=R.subset(tw,tw+dT)*exp(MP*tw)/ZP/ZD(ilsSpe,ic,is)*M*MP*4;
		      for(size_t t=0;t<A.size();t++)
			 A[t]*=exp(M*t);
		      //A[t]/=N[t];
		      A.ave_err().write(plotsPath+"/"+name+".xmg");
		      
		      return A;
		    };
		  
		  const djvec_t cV0=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_V0_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDS);
		  const djvec_t cdthV=
		    insOnDs?
		    g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_V1_TW"+to_string(itw)+"_mc"+to_string(ic)+"_dth",0,ZDS):
		    g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_dth_V1_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDS);
		  const djvec_t cdthV_unsme=
		    insOnDs?
		    g("S_m"+lsSpe+"__S_m"+lsSeq+"_V1_TW"+to_string(itw)+"_mc"+to_string(ic)+"_dth",0,ZDL):
		    g("S_m"+lsSpe+"__S_m"+lsSeq+"_dth_V1_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDL);
		  const djvec_t cS0=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDS);
		  
		  cdthV.ave_err().write(plotsPath+"/cdthV.xmg");
		  forward_derivative(log(cdthV)).ave_err().write(plotsPath+"/der_log_cdthV.xmg");
		  
		  grace_file_t megaeff("plots/megaeff.xmg");
		  grace_file_t megacoup("plots/megacoup.xmg");
		  megaeff.set_color_scheme({grace::RED,grace::BLUE,grace::BLACK,grace::GREEN4});
		  for(const auto& [tag1,tag2] : {std::pair{"V0",""},{"S0",""},{insOnDs?"V1":"dth_V1",insOnDs?"_dth":""}})
		    for(const string& t : {"_sme"+to_string(is),string("")})
		      {
			const djvec_t v=getP5P5("S_m"+lsSpe+t+
						"__S_m"+lsSeq+"_"+tag1+"_TW"+to_string(itw)+"_mc"+to_string(ic)+tag2,0);
			
			djvec_t e(v.size());
			for(size_t t=0;t<e.size();t++)
			  e[t]=log(v[(t+tw)%T]/v[(t+1+tw)%T]);
			
			djvec_t c(e.size());
			for(size_t t=0;t<e.size();t++)
			  c[t]=v[(t+tw)%T]/exp(-e[t]*t);
			
			if(t!="")
			  c*=ZDL(0,ic,0).ave()/ZDS(0,ic,0).ave();
			
			megaeff.write_vec_ave_err(e.ave_err());
			megaeff.set_legend(tag1+(std::string)tag2+t);
			
			megacoup.write_vec_ave_err(c.ave_err());
			megacoup.set_legend(tag1+(std::string)tag2+t);
		      }
		  const djvec_t D_sme_sme=getP5P5("S_m"+lsSpe+"_sme"+to_string(is)+"__sme"+to_string(is)+"_S_mc"+to_string(ic),0).symmetrized();
		  megaeff.write_vec_ave_err(effective_mass(D_sme_sme).ave_err());
		  
		  const string suff="_"+Dname+"_sme"+to_string(is)+"_tw"+to_string(tw)+"to"+Pname+".xmg";
		  
		  // const djack_t cc=ZDS(ilsSpe,ic,is)/ZDL(ilsSpe,ic,is);
		  const djvec_t cV0Unsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_V0_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDL);
		  const djvec_t cS0Unsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0,ZDL);
		  const djack_t dM=
		    //constant_fit(-forward_derivative(log(cdthV-cdthVUnsme)),18,21,"plots/dsme_nsme_S_m"+lsSpe+"__S_m"+lsSeq+"_dth_V1_TW"+to_string(itw)+"_mc"+to_string(ic)+".xmg");
		    constant_fit(-forward_derivative(log(cdthV-cdthV_unsme)),18,21,"plots/dsme_nsme_S_m"+lsSpe+"__S_m"+lsSeq+"_V1_TW"+to_string(itw)+"_mc"+to_string(ic)+"_dth.xmg");
		  
		  // const djvec_t cdthTUnsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_dth_T1_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  // const djvec_t cS0Unsme=g("S_m"+lsSpe+"__S_m"+lsSeq+"_S0_TW"+to_string(itw)+"_mc"+to_string(ic),0)*cc;
		  
		  const djvec_t& cf0=cf0List(itw,ilsCombo,ic,isme)=cS0*(mc-mseq)/(sqr(M)-sqr(MP));
		  cout<<cf0.size()<<endl;
		  f0[ic]=constant_fit(cf0,tmin3,tmax3);
		  forward_derivative(-log(cf0)).ave_err().write(plotsPath+"/effMassF0"+suff);
		  
		  const djvec_t cZ=cS0/cV0*(mc-mseq)/(M-MP);//*0+0.7;
		  cZ.ave_err().write(plotsPath+"/Z"+suff);
		  // const djack_t ZRen=cZ[30];
		  
		  const djvec_t cfp_sme=cfpList(itw,ilsCombo,ic,isme)=
		    insOnDs?
		    ((cV0+(MP-M)*cdthV)*cZ/(2*MP)):
		    ((cV0+(M-MP)*cdthV)*cZ/(2*M));
		  constant_fit(cfp_sme,tmin3,tmax3,plotsPath+"/fp_sme"+suff);
		  const djvec_t& cfp=cfpList(itw,ilsCombo,ic,isme)=
		    insOnDs?
		    ((cV0+(MP-M)*cdthV_unsme)*cZ/(2*MP)):
		    ((cV0+(M-MP)*cdthV_unsme)*cZ/(2*M));
		  fp[ic]=constant_fit(cfp,tmin3,tmax3,plotsPath+"/fp"+suff);
		  
		  const djvec_t cfpBare=
		    insOnDs?
		    ((cV0+(MP-M)*cdthV)/(2*MP)):
		    ((cV0+(M-MP)*cdthV)/(2*M));
		  cfpBare.ave_err().write(plotsPath+"/bare_fp"+suff);
		  forward_derivative(cfpBare).ave_err().write(plotsPath+"/der_bare_fp"+suff);
		  
		  const djvec_t cfpBareUnsme=
		    insOnDs?
		    ((cV0+(MP-M)*cdthV_unsme)/(2*MP)):
		    ((cV0+(M-MP)*cdthV_unsme)/(2*M));
		  cfpBareUnsme.ave_err().write(plotsPath+"/bare_fp_unsme"+suff);
		  
		  // const djvec_t cfpUnsme=(cV0Unsme+(M-MP)*cdthVUnsme)*ZRen/(2*M);
		  
		  // auto yy=[&getP5P5,is,lsSeq,itw,ic,lsSpe,suff,&Dname](const string& tag1,
		  // 						       const string& tag2="")
		  // {
		  //   const djvec_t& R=getP5P5("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_"+tag1+"_TW"+to_string(itw)+"_mc"+to_string(ic)+tag2,0);
		  //   djack_t Z,M;
		  //   two_pts_fit(Z,M,R,T/2,45,60,plotsPath+"/fit_"+tag1+tag2+"_"+Dname+"_sme"+to_string(is)+".xmg","",0);
		    
		  //   djvec_t S(R.size());
		  //   for(size_t t=0;t<T;t++)
		  //     S[t]=R[t]-two_pts_corr_fun(Z,M,T/2.0,t,0);
		    
		  //   (effective_mass(S,T/2,0)-M).ave_err().write(plotsPath+"/excFit_"+tag1+tag2+"_"+Dname+"_sme"+to_string(is)+".xmg");
		  // };
		  
		  // yy("V0");
		  // yy("V1","_dth");
		  // grace_file_t fpCompa(plotsPath+"/fpCompa"+suff);
		  // fpCompa.write_vec_ave_err(cfp.ave_err());
		  // fpCompa.set_legend("Smeared");
		  // fpCompa.set_all_colors(grace::BLUE);
		  // fpCompa.write_vec_ave_err(cfpUnsme.ave_err());
		  // fpCompa.set_legend("Local");
		  
		  jack_fit_t fitter;
		  djvec_t pars(3);
		  fitter.add_fit_par(pars[0],"C",4,0.1);
		  fitter.add_fit_par(pars[1],"dM",MPList[0].ave(),0.1);
		  fitter.add_fit_par(pars[2],"dZ",-2.6,0.1);
		  
		  const size_t t2s=11;
		  for(size_t it=t2s;it<=dT;it++)
		    fitter.add_point(cfp[it],
				     [t=(double)it](const vector<double>& pars,
						    const size_t iel)
				     {
				       return pars[0]+pars[2]*exp(-pars[1]*t);
				     });
		  fitter.fit();
		  for(const auto& p : pars)
		    cout<<""<<p.ave_err()<<endl;
		  
		  cout<<"dM: "<<(pars[1]/a).ave_err()<<endl;
		  
		  grace_file_t plot2sFitFp(plotsPath+"/2sfp"+suff);
		  plot2sFitFp.write_vec_ave_err(cfp.ave_err());
		  for(const auto& [b,e] : {std::make_pair(0,t2s),{t2s,dT}})
		    plot2sFitFp.write_polygon([&pars](const double& t)
		    {
		      return pars[0]+pars[2]*exp(-pars[1]*t);
		    },b,e);
		  
		  // grace_file_t plot2sFit2subFp(plotsPath+"/2sfp2sub"+suff);
		  // djvec_t i=cfp;
		  // for(size_t t=0;t<cfp.size();t++)
		  //   i[t]-=pars[2]*exp(-pars[1]*t);
		  // plot2sFit2subFp.write_vec_ave_err(i.ave_err());
		  // plot2sFit2subFp.write_polygon([&pars](const double& t)
		  //   {
		  //     return pars[0];
		  //   },0,dT);
		  
		  // cout<<" t dom: "<<(-log(abs(pars[0]/pars[2]/100))/pars[1]).ave_err()<<endl;
		  
		  // const djvec_t cfpfr0=(cV0+(M-MP)*cdthV)*cZ/(2*M)/cf0;
		  // constant_fit(cfpfr0,tmin3,tmax3,plotsPath+"/fpfr0"+suff);
		  
		  djvec_t cft;
		  
		  // const djvec_t cdthT=g("S_m"+lsSpe+"_sme"+to_string(is)+"__S_m"+lsSeq+"_dth_T1_TW"+to_string(itw)+"_mc"+to_string(ic),0);
		  // cft=cdthT*(sqr(M)+sqr(MP))/(2*M);
		  // ft[ic]=constant_fit(cft,tmin3,tmax3,plotsPath+"/ft"+suff);
		  
		  for(const auto&[tag,c,f] :
			{make_tuple("p",cfp,fp[ic]),
			 /*      */{"0",cf0,f0[ic]},
			 /*      */{"t",cft,ft[ic]}
			})
		    {
		      grace_file_t plotf(plotsPath+"/f"+(tag+suff));
		      plotf.set_title((string)tag+"("+Dname+"->"+Pname+")");
		      plotf.set_xaxis_label("t/a");
		      plotf.write_vec_ave_err(c.ave_err());
		      plotf.write_constant_band(tmin3,tmax3,f);
		      for(const string coll : {"etm","fnal"})
			if(const auto i=liter.find(coll+tag+Dname+"To"+Pname);i!=liter.end())
			  {
			    plotf.set_legend(coll);
			    plotf.write_constant_band(0,dT,djack_t(gauss_filler_t(i->second,12432)));
			  }
			// else
			//   cout<<coll+tag+Dname+"To"+Pname<<" not found in the literature"<<endl;
		    }
		  
		  cfp.bin_write("cfp.dat");
		  cf0.bin_write("cf0.dat");
		}
	    }
	}
      
      const string Dname=
	(string)"D"+((ilsSpe==1)?"s":"");
      auto plot=
	[&](const string& tag,
	      const djvec_t& f)
	{
	  const string decay=tag+Dname+"to"+Pname;
	  grace_file_t plotf(plotsPath+"/f"+decay+".xmg");
	  for(size_t ic=0;ic<nmc;ic++)
	    plotf.write_ave_err(MD(0,ic,0).ave(),f[ic].ave_err());
	};
      
      for(const auto&[tag,f] : {make_pair("p",fp),{"0",f0},{"t",ft}})
	plot(tag,f);
    }
  
  // {
  //   djvec_t& f=cfpList(0,0,0,0);
    
  //   jack_fit_t fitter;
  //   djvec_t pars(3);
  //   const double guessDelta=delta(0,0,0,0)[20].ave();
  //   fitter.add_fit_par(pars[0],"A",f[20].ave(),0.1);
  //   fitter.add_fit_par(pars[1],"B",-0.1,0.1);
  //   fitter.add_fit_par(pars[2],"C",guessDelta,0.1);
    
  //   cout<<"Fitting exponential on top of fp"<<endl;
  //   for(size_t it=tpExpMin;it<=tpExpMax;it+=2)
  //     fitter.add_point([num=f[it]](const vector<double> &p,int iel){return num[iel];},
  // 		       [t=(double)it](const vector<double>& pars,
  // 				      const size_t iel)
  // 		       {
  // 			 return pars[0]+pars[1]*exp(-pars[2]*t);
  // 		       },f[it],0);
  //   const auto [finalCh2,nDof,status]=fitter.fit(true);
    
  //   for(const auto& p : pars)
  //     cout<<" "<<p.ave_err()<<endl;
  //   cout<<finalCh2.ave()<<endl;
    
  //   cout<<"Mfit: "<<pars[2].ave_err()<<" significativity of difference with delta "<<guessDelta<<": "<<(pars[2]-guessDelta).significativity()<<endl;
    
  //   effective_mass(f-pars[0],T/2,0).ave_err().write(plotsPath+"/fpExc.xmg");
    
  //   grace_file_t plotFpDiffFit(plotsPath+"/fpDiff.xmg");
  //   plotFpDiffFit.write_vec_ave_err(f.ave_err());
  //   plotFpDiffFit.write_polygon([&pars](const double& t)
  //   {
  //     return pars[0]+pars[1]*exp(-pars[2]*t);
  //   },tpExpMin,tpExpMax);
  // }
  
  for(auto& [tag,c] :
	{make_tuple("p",&cfpList),
	 {"0",&cf0List}})
    for(size_t ic=0;ic<nmc;ic++)
      for(size_t ilsCombo=0;ilsCombo<lsCombo.size();ilsCombo++)
	{
	  const size_t& ilsSeq=lsCombo[ilsCombo].first;
	  const size_t& ilsSpe=lsCombo[ilsCombo].second;
	  
	  const string Pname=
	    (ilsSpe==1 or ilsSeq==1)?"K":"Pi";
	  const string Dname=
	    "D"+to_string(ic)+((ilsSpe==1)?"s":"");
	  const string suff="_"+Dname+"_sme"+to_string(isme)+"to"+Pname+".xmg";
	  
	  vector<double> x(dT+1);
	  for(size_t t=0;t<=dT;t++)
	    x[t]=t*a.ave();
	  
	  grace_file_t plotf(plotsPath+"/f"+(tag+suff));
	  plotf.set_title((string)tag+"("+Dname+"->"+Pname+")");
	  plotf.set_xaxis_label("t/a");
	  
	  const std::vector<grace::color_t> colors{grace::BLACK,grace::RED,grace::BLUE};
	  for(size_t itw=0;itw<ntw;itw++)
	    {
	      const size_t& tw=twList[itw];
	      const vec_ave_err_t y=(*c)(itw,ilsCombo,ic,isme).ave_err();
              plotf.write_vec_ave_err(x,y);
	      plotf.set_legend(std::to_string(tw));
	      plotf.set_all_colors(colors[itw%colors.size()]);
	    }
	  
	  // const djvec_t cfpExtr=
	  //   2*(*c)(1,ilsCombo,ic,isme)-(*c)(0,ilsCombo,ic,isme);
	  
	  djvec_t cfExtr(dT);
	  
	  vector<double> w;
	  for(size_t itw=0;itw<ntw;itw++)
	    w.push_back(exp(2*MPList[0].ave()*twList[itw]));
	  const double wmax=
	    *std::max_element(w.begin(),w.end());
	  
	  if(ntw>1)
	    {
	      djvec_t linExtrTerm(dT);
	      for(size_t t=0;t<dT;t++)
		{
		  djvec_t y(ntw);
		  for(size_t itw=0;itw<ntw;itw++)
		    y[itw]=(*c)(itw,ilsCombo,ic,isme)[t];
		  
		  const djvec_t p=poly_fit(w,y,1);
		  // cout<<chi2_poly_fit(x,y,1,0.0,xmax,p)<<endl;
		  write_poly_fit_plot(plotsPath+"/extr"+(string)tag+"_"+to_string(t)+".xmg",0.0,wmax,p,w,y);
		  cfExtr[t]=p[0];
		  linExtrTerm[t]=p[1];
		  
		  // cout<<t<<" "<<(p[1]/p[0]/exp(-MPList[0].ave()*T)).ave_err()<<endl;
		}
	      plotf.write_vec_ave_err(x,cfExtr.ave_err());
	      plotf.set_legend("extr");
	      plotf.set_all_colors(grace::GREEN4);
	      
	      const size_t tExcFitMin=15;
	      const size_t tExcFitMax=25;
	      const djack_t dM=constant_fit(symmetric_derivative(-log(symmetric_derivative(cfExtr))),tExcFitMin,tExcFitMax,(string)plotsPath+"/derLogDerF"+tag+suff);
	      djvec_t dcF=-log(symmetric_derivative(cfExtr));
	      for(size_t t=0;t<dT;t++)
		dcF[t]-=dM*t;
	      const djack_t dcTemp= // -log(-dc*dm))
		constant_fit(dcF,tExcFitMin,tExcFitMax,(string)plotsPath+"/excEffCouplF"+tag+suff);
	      const djack_t dc=
		-exp(-dcTemp)/dM;
	      
	      djvec_t cfSub=cfExtr;
	      for(size_t t=0;t<dT;t++)
		cfSub[t]-=dc*exp(-dM*t);
	      
	      const djack_t fSub=
		constant_fit(cfSub,tExcFitMin,tExcFitMax,(string)plotsPath+"/subF"+tag+suff);
	      
	      plotf.write_vec_ave_err(x,cfSub.ave_err());
	      plotf.set_legend("sub");
	      plotf.set_all_colors(grace::VIOLET);
	      
	      for(const string coll : {"etm","fnal"})
		if(const auto i=liter.find(coll+tag+Dname+"To"+Pname);i!=liter.end())
		  {
		    plotf.set_legend(coll);
		    plotf.write_constant_band(0,dT*a.ave(),djack_t(gauss_filler_t(i->second,12432)));
		  }
	      plotf.new_data_set();
	      
	      plotf.set_legend("extrSubFit");
	      plotf.write_polygon([&fSub,
				   &dc,
				   &dM,
				   aave=a.ave()](const double& t)
	      {
		return fSub+dc*exp(-dM*t/aave);
	      },
				  1,dT*a.ave(),grace::ORANGE);
	      
	      linExtrTerm.ave_err().write((string)plotsPath+"/f"+tag+"LinTerm.xmg");
	      symmetric_derivative(linExtrTerm).ave_err().write((string)plotsPath+"/f"+tag+"LinTermDer.xmg");
	      symmetric_derivative(log(symmetric_derivative(linExtrTerm))).ave_err().write((string)plotsPath+"/f"+tag+"LinTermDerLogDer.xmg");
	    }
	}
  
  return 0;
}
