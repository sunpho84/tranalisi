#include <set>
#include <tranalisi.hpp>

#include <unsupported/Eigen/MatrixFunctions>

const size_t T=128;
const size_t nCombo=225;
const size_t nG=4;

//! perform a fit to determine the slope
void two_pts_LL_SL_SS_fit(djack_t& Z_L,
			  djack_t& Z_S,
			  djack_t& M,
			  const djvec_t &corr_LL,
			  const djvec_t &corr_SL,
			  const djvec_t &corr_SS,
			  const size_t tmin,
			  const size_t tmax,
			  const string& path)
{
  const size_t TH=T/2;
  
  //perform a preliminary fit
  djack_t Z_LL,Z_SS;
  two_pts_fit(Z_SS,M,corr_SS,TH,tmin,tmax,"/tmp/test_eff_mass_sme_sme.xmg");
  two_pts_fit(Z_LL,M,corr_LL,TH,tmin,tmax,"/tmp/test_eff_mass_loc_loc.xmg");
  Z_L=sqrt(Z_LL);
  Z_S=sqrt(Z_SS);
  
  //parameters to fit
  minimizer_pars_t pars;
  pars.add("M",M[0],M.err());
  pars.add("Z_L",Z_L[0],Z_L.err());
  pars.add("Z_S",Z_S[0],Z_S.err());
  
  //! fit for real
  size_t iel=0;
  auto x=vector_up_to<double>(corr_SL.size());
  multi_ch2_t<djvec_t> two_pts_SL_fit_obj({x,x,x},{tmin,tmin,tmin},{tmax,tmax,tmax},{corr_LL,corr_SL,corr_SS},
					  {two_pts_corr_fun_t(TH,1),
					   two_pts_corr_fun_t(TH,1),
					   two_pts_corr_fun_t(TH,1)},
					  [](const vector<double> &p,size_t icontr)
					  {
					    switch(icontr)
					      {
					      case 0:return vector<double>({p[1]*p[1],p[0]});break;
					      case 1:return vector<double>({p[2]*p[1],p[0]});break;
					      case 2:return vector<double>({p[2]*p[2],p[0]});break;
					      default: CRASH("Unknown contr %zu",icontr);return p;
					      }
					  },iel);
      
  minimizer_t minimizer(two_pts_SL_fit_obj,pars);
  
  for(iel=0;iel<corr_SL[0].size();iel++)
    {
      //minimize and print the result
      vector<double> par_min=minimizer.minimize();
      M[iel]=par_min[0];
      Z_L[iel]=par_min[1];
      Z_S[iel]=par_min[2];
    }
  
  //write plots
  if(path!="")
    {
      grace_file_t out(path);
      out.write_polygon([&M](double x){return M;},tmin,tmax);
      out.write_vec_ave_err(x,effective_mass(corr_LL,TH).ave_err());
      out.write_vec_ave_err(x,effective_mass(corr_SL,TH).ave_err());
      out.write_vec_ave_err(x,effective_mass(corr_SS,TH).ave_err());
    }
}

int main()
{
  set_njacks(15);
  
  raw_file_t fin("data","r");
  
  vector<double> raw;
  
  vector<array<size_t,4>> comboTag(nCombo);
  
  vector<string> rawComboTag;
  
  const string compressed_data_path="compressed_data";
  if(not file_exists(compressed_data_path))
    {
      bool comboRead=false;
      while(not fin.feof())
	{
	  for(size_t iC=0;iC<nCombo;iC++)
	    {
	      string A,B;
	      fin.expect("#");
	      fin.expect("Contraction");
	      fin.expect("of");
	      fin.read(A);
	      //cout<<A<<endl;
	      fin.expect("^");
	      fin.expect("\\dag");
	      fin.expect("and");
	      fin.read(B);
	      //cout<<B<<endl;
	      
	      if(not comboRead)
		{
		  rawComboTag.push_back(A+" "+B);
		  auto dec=
		    [](const string& s) -> pair<size_t,size_t>
		    {
		      size_t a{},b{};
		      sscanf(s.c_str(),"%zu_S_%zu",&a,&b);
		      
		      return {a,b};
		    };
		  
		  const auto [A1,A2]=dec(A);
		  const auto [B1,B2]=dec(B);
		  
		  comboTag[iC]={A1,A2,B1,B2};
		}
	      
	      for(size_t iG=0;iG<nG;iG++)
		{
		  string C;
		  fin.expect("#");
		  fin.read(C);
		  //cout<<C<<endl;
		  
		  for(size_t t=0;t<T;t++)
		    raw.push_back(fin.read<double>());
		}
	      
	      fin.skip_line();
	      fin.skip_line();
	      fin.skip_line();
	    }
	  
	  comboRead=true;
	}
      
      raw_file_t fout(compressed_data_path,"w");
      fout.bin_write(raw.size());
      fout.bin_write(raw);
      fout.bin_write(comboTag.size());
      fout.bin_write(comboTag);
      fout.bin_write(rawComboTag.size());
      fout.bin_write(rawComboTag);
    }
  else
    {
      raw_file_t fin(compressed_data_path,"r");
      raw.resize(fin.bin_read<size_t>());
      fin.bin_read(raw);
      comboTag.resize(fin.bin_read<size_t>());
      fin.bin_read(comboTag);
      rawComboTag.resize(fin.bin_read<size_t>());
      fin.bin_read(rawComboTag);
    }
  
  const size_t nConfs=raw.size()/T/4/nCombo;
  cout<<"NConfs: "<<nConfs<<endl;
  
  index_t idRaw({{"Conf",nConfs},{"Combo",nCombo},{"Gamma",nG},{"T",T}});
  index_t id({{"Conf",nConfs},{"Combo",nCombo},{"T",T/2+1}});
  
  vector<double> BKBK(id.max());
  for(size_t iConf=0;iConf<nConfs;iConf++)
    for(size_t iC=0;iC<nCombo;iC++)
      for(size_t t=0;t<T/2+1;t++)
	{
	  double& d=BKBK[id({iConf,iC,t})];
	  d=0;
	  for(size_t iG=1;iG<nG;iG++)
	    d+=(raw[idRaw({iConf,iC,iG,t})]+raw[idRaw({iConf,iC,iG,(T-t)%T})])/2;
	  d/=3;
	  
	  //d=(raw[idRaw({iConf,iC,0,t})]+raw[idRaw({iConf,iC,0,(T-t)%T})])/2;
	}
  
  /// Creates the single smearings
  map<size_t,size_t> singleSme;
  for(const auto& [A1,A2,A3,A4] : comboTag)
    {
      singleSme[A1+A3];
      singleSme[A2+A4];
    }
  
  /// Number of single smearing
  const size_t nSme=singleSme.size();
  
  /// Label smearing
  size_t iSingle=0;
  for(auto& [A,B] : singleSme)
    {
      B=iSingle++;
      cout<<A<<" "<<B<<endl;
    }
  
  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
  
  const size_t INVALID=999999999;
  vector<size_t> comboMap(nCombo,INVALID);
  
  for(size_t iC=0;iC<nCombo;iC++)
    {
      const auto& [A1,A2,A3,A4]=comboTag[iC];
      
      const size_t A=singleSme[A1+A3];
      const size_t B=singleSme[A2+A4];
      
      //cout<<"combo "<<iC<<" "<<rawComboTag[iC]<<" corresponds to "<<A2+A3<<" "<<A1+A4<<" translated into "<<A<<" "<<B<<endl;
      
      comboMap[B+nSme*A]=iC;
      
    }
  
  grace_file_t effPlot("eff.xmg");
  vector<djvec_t> C(nCombo,djvec_t(T/2+1,0));
  for(size_t iC=0;iC<nCombo;iC++)
    {
      jackknivesFill(nConfs,
		     [&C,
		      &iC,
		      &id,
		      &BKBK](const size_t& iConf,
			     const size_t& iClust,
			     const double& w)
		     {
		       for(size_t t=0;t<T/2+1;t++)
			 C[iC][t][iClust]+=BKBK[id({iConf,iC,t})]*w;
		     });
      C[iC].clusterize((double)nConfs/njacks);
    }
  
  for(size_t iSme=0;iSme<nSme;iSme++)
    {
      const size_t iC=comboMap[iSme+nSme*iSme];
      if(iC==INVALID)
	CRASH("searching for invalid combo %zu %zu",iSme,iSme);
      effPlot.write_vec_ave_err(effective_mass(C[iC]).ave_err());
      
      string legend=std::to_string(iC)+" ";
      for(const auto& i : comboTag[iC])
	legend+=std::to_string(i)+" ";
      effPlot.set_legend(legend);
      effPlot.new_data_set();
    }
  
  /////////////////////////////////////////////////////////////////
  
  //Fill the rhs matrix
  input_file_t input("input.txt");
  [[maybe_unused]] const size_t t0=input.read<size_t>("t0");
  const size_t nLev=input.read<size_t>("NLev");
  vector<size_t> smeToUse(nLev);
  for(size_t iLev=0;iLev<nLev;iLev++)
    {
      input.read(smeToUse[iLev]);
      if(singleSme.find(smeToUse[iLev])==singleSme.end())
	CRASH("Unable to find sme %zu",smeToUse[iLev]);
    }
  
#if 1
  typedef Matrix<double,Dynamic,Dynamic> Matr;
  
  vector<djvec_t> cForGevp(nLev*nLev);
  for(size_t iLev=0;iLev<nLev;iLev++)
    for(size_t jLev=0;jLev<nLev;jLev++)
      {
	const auto [niSme,njSme]=minmax(smeToUse[iLev],smeToUse[jLev]);
	const size_t iSme=singleSme[niSme];
	const size_t jSme=singleSme[njSme];
	const size_t iC=comboMap[jSme+nSme*iSme];
	
	cout<<"for "<<iLev<<" "<<jLev<<" corresponding to "<<niSme<<" "<<njSme<<" pointing to "<<iSme<<" "<<jSme<<" using "<<iC<<" "<<rawComboTag[iC]<<endl;
	if(iC==INVALID)
	  CRASH("invalid");
	
	cForGevp[iLev*nLev+jLev]=C[iC];
      }
  
  GeneralizedEigenSolver<Matr> ge;
  vector<djvec_t> eig(nLev,djvec_t(T/2+1));
  vector<djvec_t> w(nLev*nLev,djvec_t(T/2+1));
  for(size_t ijack=0;ijack<=njacks;ijack++)
    {
      //Fill the rhs matrix
      Matr b0(nLev,nLev);
      for(size_t i=0;i<nLev;i++)
	for(size_t j=0;j<nLev;j++)
	  b0(i,j)=cForGevp[i*nLev+j][t0][ijack];
      
      for(size_t t=0;t<T/2+1;t++)
	{
	  Matr b(nLev,nLev);
	  for(size_t i=0;i<nLev;i++)
	    for(size_t j=0;j<nLev;j++)
	      b(i,j)=cForGevp[i*nLev+j][t][ijack];
	  
	  ge.compute(b,b0);
	  
	  for(size_t iLev=0;iLev<nLev;iLev++)
	    eig[iLev][t][ijack]=ge.eigenvalues()(iLev).real();
	  
	  // const Matr o=(ge.eigenvectors().adjoint()*b*ge.eigenvectors()).real();
	  // cout<<o<<endl;
	  
	  for(size_t iLev=0;iLev<nLev;iLev++)
	    for(size_t iOp=0;iOp<nLev;iOp++)
	      w[iOp+nLev*iLev][t][ijack]=ge.eigenvectors()(iOp,iLev).real();
	}
    }
  
  vector<vector<size_t>> order(T/2+1,vector<size_t>(nLev));
  for(size_t t=0;t<T/2+1;t++)
    {
      vector<pair<double,size_t>> d(nLev);
      
      for(size_t iLev=0;iLev<nLev;iLev++)
	d[iLev]={eig[iLev][t][njacks],iLev};
      
      sort(d.begin(),d.end(),std::greater{});
      
      for(size_t iLev=0;iLev<nLev;iLev++)
	order[t][iLev]=d[iLev].second;
    }
  
  grace_file_t wFile("w.xmg");
  const size_t lev0=0;
  for(size_t iOp=0;iOp<nLev;iOp++)
    {
      djvec_t y(T/2+1);
      for(size_t t=0;t<T/2+1;t++)
	{
	  const double s=(w[0+nLev*order[t][lev0]][t][0]>0)?-1:+1;
	  y[t]=w[iOp+nLev*order[t][lev0]][t]*s;
	}
      
      wFile.write_vec_ave_err(y.ave_err());
    }
  
  grace_file_t gaussFile("gauss.xmg");
  const size_t t1=10;
  const double k=0.4;
  const double xmax=9.0;
  for(size_t iOp=0;iOp<nLev;iOp++)
    gaussFile.write_line([sig=sqrt(smeToUse[iOp]*k/(1+6*k))](const double& x)
    {
      return gauss(0.0,sig,x);
    },0,xmax);
  gaussFile.write_line([&smeToUse,
			&k,
			&nLev,
			&w,
			&order](const double& x)
  {
    double s=0;
    for(size_t iOp=0;iOp<nLev;iOp++)
      s+=gauss(0.0,sqrt(smeToUse[iOp]*k/(1+6*k)),x)*w[iOp+nLev*order[t1][lev0]][t1][njacks];
    return s;
  },
		       0,xmax);
  
  cout<<"/////////////////////////////////////////////////////////////////"<<endl;
  
  grace_file_t eigFile("eig.xmg");
  for(size_t iLev=0;iLev<nLev;iLev++)
    eigFile.write_vec_ave_err(effective_mass(eig[iLev]).ave_err());
  
  djvec_t reco(T/2+1);
  reco=0.0;
  for(size_t t=0;t<T/2+1;t++)
    for(size_t iOp=0;iOp<nLev;iOp++)
      for(size_t jOp=0;jOp<nLev;jOp++)
	reco[t]+=w[iOp+nLev*order[t][lev0]][t1]*cForGevp[iOp*nLev+jOp][t]*w[jOp+nLev*order[t][lev0]][t1];
  effective_mass(reco).ave_err().write("reco.xmg");
  
  grace_file_t oriFile("ori.xmg");
  for(size_t iLev=0;iLev<nLev;iLev++)
    for(size_t jLev=iLev;jLev<nLev;jLev++)
      {
	const auto [niSme,njSme]=minmax(smeToUse[iLev],smeToUse[jLev]);
	const size_t iSme=singleSme[niSme];
	const size_t jSme=singleSme[njSme];
	const size_t iC=comboMap[jSme+nSme*iSme];
	
	cout<<iC<<endl;
	
	oriFile.write_vec_ave_err(effective_mass(C[iC]).ave_err());
	oriFile.set_legend(rawComboTag[iC]);
	oriFile.need_close_set=true;
      }
#endif
  
  /////////////////////////////////////////////////////////////////
  
  grace_file_t checkFile("check.xmg");
  for(size_t iLev=0;iLev<nLev;iLev++)
    for(size_t jLev=iLev+1;jLev<nLev;jLev++)
      {
	const auto [niSme,njSme]=minmax(smeToUse[iLev],smeToUse[jLev]);
	const size_t iSme=singleSme[niSme];
	const size_t jSme=singleSme[njSme];
	const size_t iC11=comboMap[iSme+nSme*iSme];
	const size_t iC12=comboMap[jSme+nSme*iSme];
	const size_t iC22=comboMap[jSme+nSme*jSme];
	
	cout<<"Combining C["<<iC12<<"]/sqrt(C["<<iC11<<"]*C["<<iC22<<"])"<<endl;
	
	checkFile.write_vec_ave_err((C[iC12]/sqrt(C[iC11]*C[iC22])).ave_err());
	
	djvec_t A11=C[iC11],A12=C[iC12],A22=C[iC22];
	
	// djack_t Z_L,Z_S,M;
	// two_pts_LL_SL_SS_fit(Z_L,Z_S,M,A11,A12,A22,35,50,"LL_SL_SS_fit.xmg");
	
	djack_t ZLL,ZLS,ZSS,M_LL,M_LS,M_SS;
	two_pts_fit(ZLL,M_LL,A11,T/2,35,50,"two_pts_LL.xmg");
	two_pts_fit(ZLS,M_LS,A12,T/2,35,50,"two_pts_LS.xmg");
	two_pts_fit(ZSS,M_SS,A22,T/2,35,50,"two_pts_SS.xmg");
	
	cout<<"sqrt(ZLL*ZSS): "<<smart_print(sqrt(ZLL*ZSS))<<"     vs    ZLS: "<<smart_print(ZLS)<<endl;
	cout<<"M_LL: "<<M_LL.ave_err()<<endl;
	cout<<"M_LS: "<<M_LS.ave_err()<<endl;
	cout<<"M_SS: "<<M_SS.ave_err()<<endl;
	
	for(size_t t=0;t<T/2+1;t++)
	  {
	    A11[t]-=two_pts_corr_fun(ZLL,M_LL,T/2.0,t,1);
	    A12[t]-=two_pts_corr_fun(ZLS,M_LS,T/2.0,t,1);
	    A22[t]-=two_pts_corr_fun(ZSS,M_SS,T/2.0,t,1);
	  }
	
	djack_t Zb_L,Zb_S,Mb;
	two_pts_LL_SL_SS_fit(Zb_L,Zb_S,Mb,A11,A12,A22,20,35,"LL_SL_SS_fit_b.xmg");
	
	cout<<iSme<<" "<<jSme<<" "<<(sqrt(ZLL*ZSS)/ZLS-1).ave_err()<<"    AAA"<<endl;
	
      }
  
  return 0;
}
