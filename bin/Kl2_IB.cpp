#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

int T;
index_t<4> ind;
size_t noa=8,nbeta=3,nboots=100;

class lat_par_t
{
public:
  dboot_t ml,ms,mc;
  dbvec_t ainv(nbeta);
};

dboot_t read_boot(const raw_file_t &file)
{
  dboot_t out;
  for(size_t ib=0;ib<nboots;ib++) file.read(out[ib]);
  return out;
}

djvec_t load(const char *path,size_t im1,size_t im2,size_t reim=RE,int rpar=1,int spat_par=1)
{
  djvec_t corr_r0=read_djvec(path,T,ind({im1,im2,0,reim}));
  djvec_t corr_r1=read_djvec(path,T,ind({im1,im2,1,reim}));

  return djvec_t(corr_r0+rpar*corr_r1).symmetrized(spat_par)*0.5;
}
  
int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s T",arg[0]);
  
  T=atoi(arg[1]);
  int TH=T/2;
  
  ind.set_ranges({3,3,2,2});  
  set_njacks(15);
  
  //pion
  
  djvec_t pi_ratio=load("corrLL_P5P5",0,0)/load("corr00_P5P5",0,0);

  pi_ratio.ave_err().write("pi_ratio.xmg");

  djack_t pi_M,pi_A,pi_SL;
  two_pts_with_ins_ratio_fit(pi_M,pi_A,pi_SL,load("corr00_P5P5",0,0),load("corrLL_P5P5",0,0),TH,10,20,"Mpi.xmg","pi_ins.xmg");
  cout<<"pi_M: "<<pi_M.ave_err()<<endl;
  cout<<"pi_A: "<<pi_A.ave_err()<<endl;
  cout<<"pi_SL: "<<pi_SL.ave_err()<<endl;


  raw_file_t obs_file("obs_pi","w");
  obs_file.bin_write(pi_M);
  obs_file.bin_write(pi_SL);

  
  //kaon


  djvec_t k_ratio_exch=load("corrLL_P5P5",0,1)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_self=load("corr0M_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_tad=load("corr0T_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_s=load("corr0S_P5P5",1,0)/load("corr00_P5P5",1,0);
  djvec_t k_ratio_p=load("corr0P_P5P5",1,0,IM,-1)/load("corr00_P5P5",1,0);

  k_ratio_exch.ave_err().write("k_ratio_exch.xmg");
  k_ratio_self.ave_err().write("k_ratio_self.xmg");
  k_ratio_tad.ave_err().write("k_ratio_tad.xmg");
  k_ratio_s.ave_err().write("k_ratio_s.xmg");
  k_ratio_p.ave_err().write("k_ratio_p.xmg");



  djack_t k_M1,k_M2,k_M3,k_M4,k_A_exch,k_SL_exch,k_A_selftad,k_SL_selftad,k_A_s,k_SL_s,k_A_p,k_SL_p;

  djack_t k_M,k_Z;

  djvec_t corr_00=load("corr00_P5P5",1,0);
  two_pts_migrad_fit(k_Z,k_M,corr_00,TH,10,20,"kaon_mass.xmg");
  cout<<"k_Z: "<<k_Z.ave_err()<<endl;
  cout<<"k_M: "<<k_M.ave_err()<<endl;

  
  
  two_pts_with_ins_ratio_fit(k_M1,k_A_exch,k_SL_exch,load("corr00_P5P5",1,0),load("corrLL_P5P5",0,1),TH,10,20,"Mk1.xmg","k_exch.xmg");
  two_pts_with_ins_ratio_fit(k_M2,k_A_selftad,k_SL_selftad,load("corr00_P5P5",1,0),djvec_t(load("corr0M_P5P5",1,0)+load("corr0T_P5P5",1,0)),TH,10,20,"Mk2.xmg","k_selftad.xmg");
  two_pts_with_ins_ratio_fit(k_M3,k_A_s,k_SL_s,load("corr00_P5P5",1,0),load("corr0S_P5P5",1,0),TH,10,20,"Mk3.xmg","k_s.xmg");
  two_pts_with_ins_ratio_fit(k_M4,k_A_p,k_SL_p,load("corr00_P5P5",1,0),load("corr0P_P5P5",1,0),TH,10,20,"Mk4.xmg","k_p.xmg");
  
  
  cout<<"k_M1: "<<k_M1.ave_err()<<endl;
  cout<<"k_A_exch: "<<k_A_exch.ave_err()<<endl;
  cout<<"k_SL_exch: "<<k_SL_exch.ave_err()<<endl;
  cout<<"k_M2: "<<k_M2.ave_err()<<endl;
  cout<<"k_A_selftad: "<<k_A_selftad.ave_err()<<endl;
  cout<<"k_SL_selftad: "<<k_SL_selftad.ave_err()<<endl;
  cout<<"k_M3: "<<k_M3.ave_err()<<endl;
  cout<<"k_A_s: "<<k_A_s.ave_err()<<endl;
  cout<<"k_SL_s: "<<k_SL_s.ave_err()<<endl;
  cout<<"k_M4: "<<k_M4.ave_err()<<endl;
  cout<<"k_A_p: "<<k_A_p.ave_err()<<endl;
  cout<<"k_SL_p: "<<k_SL_p.ave_err()<<endl;

  //D meson

 
  djvec_t D_ratio_exch=load("corrLL_P5P5",0,2)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_self=load("corr0M_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_tad=load("corr0T_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_s=load("corr0S_P5P5",2,0)/load("corr00_P5P5",2,0);
  djvec_t D_ratio_p=load("corr0P_P5P5",2,0,IM,-1)/load("corr00_P5P5",2,0);

  D_ratio_exch.ave_err().write("D_ratio_exch.xmg");
  D_ratio_self.ave_err().write("D_ratio_self.xmg");
  D_ratio_tad.ave_err().write("D_ratio_tad.xmg");
  D_ratio_s.ave_err().write("D_ratio_s.xmg");
  D_ratio_p.ave_err().write("D_ratio_p.xmg");

  raw_file_t file("/Users/mac/Programmi/Tesi_Laurea_Magistrale_bozza/Archivio definitivo/Tesi Magistrale/Programmi/New Nucleon/bootstrapFVE/ultimate_input.out","r");

  vector<lat_par_t> lat_par(noa);
  double dum;
  file.expect({"ml","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].ml=read_boot(file);
  file.expect({"ms","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].ms=read_boot(file);
  file.expect({"mc(2","GeV)","(GeV)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++) lat_par[ia].mc=read_boot(file);
  file.expect({"a^-1","(GeV)","(1.90","1.95","2.10)"});
  file.read(dum);
  for(size_t ia=0;ia<noa;ia++)
    for(size_t ibeta=0;ibeta<nbeta;ibeta++)
      for(size_t iboot=0;iboot<nboots;iboot++)
	file.read(lat_par[ia].ainv[ibeta][iboot]);
  //file.skip_line_line(2);
  //double ainv[Nbeta][Nev], r0[Nev], Z1[Nbeta], Z2[Nbeta], Zev[Nbeta][Nev], f0[Nev], B0[Nev];
  //int iboot[Nbeta][Nmusea_max][Nev];
  
  cout<<lat_par[0].ms.ave_err()<<endl;
  
  return 0;
}
