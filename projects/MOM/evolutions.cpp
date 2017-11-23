#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#define EXTERN_EVOLUTIONS
 #include <evolutions.hpp>

#define lambdaQCD 0.250
#define Z3 1.20206
#define Nc 3.0

#include <cmath>

using namespace std;

double alphas(int Nf,double ainv,double mu2)
{
  int CF;
  CF = (Nc*Nc-1.)/(2.*Nc);
  
  double lam0,L2,LL2,b1,b2,b3;
  double als0, als1, als2, als3;
  double beta_0,beta_1,beta_2,beta_3;
  
  beta_0 = (11.*Nc-2.*Nf)/3.;
  
  beta_1 = 34./3.*pow(Nc,2) - 10./3.*Nc*Nf-2*CF*Nf;
  beta_2 = (2857./54.)*pow(Nc,3) + pow(CF,2)*Nf -
    205./18.*CF*Nc*Nf -1415./54.*pow(Nc,2)*Nf +
    11./9.*CF*pow(Nf,2) + 79./54.*Nc*pow(Nf,2);
  
  beta_3 = (150653./486. - 44./9.*Z3)*pow(Nc,4) +
    (-39143./162. + 68./3.*Z3)*pow(Nc,3)*Nf +
    (7073./486. - 328./9.*Z3)*CF*pow(Nc,2)*Nf +
    (-2102./27. + 176./9.*Z3)*pow(CF,2)*Nc*Nf +
    23.*pow(CF,3)*Nf + (3965./162. + 56./9.*Z3)*pow(Nc,2)*pow(Nf,2) +
    (338./27. - 176./9.*Z3)*pow(CF,2)*pow(Nf,2) +
    (4288./243. + 112./9.*Z3)*CF*Nc*pow(Nf,2) + 53./243.*Nc*pow(Nf,3) +
    154./243.*CF*pow(Nf,3) +
    (-10./27. + 88./9.*Z3)*pow(Nc,2)*(pow(Nc,2)+36.) +
    (32./27. - 104./9.*Z3)*Nc*(pow(Nc,2)+6)*Nf +
    (-22./27. + 16./9.*Z3)*(pow(Nc,4) - 6.*pow(Nc,2) + 18.)/pow(Nc,2)*pow(Nf,2);
  
  b1=beta_1/beta_0/4./M_PI;
  b2=beta_2/beta_0/16./pow(M_PI,2);
  b3=beta_3/beta_0/64./pow(M_PI,3);
  
  lam0=lambdaQCD/ainv;
  
  L2   = log( mu2/(pow(lam0,2) ) );
  LL2  = log( L2 );
  
  als0 = 4.*M_PI/beta_0/L2;
  als1 = als0 - pow(als0,2)*b1*LL2;
  als2 = als1 + pow(als0,3)*(pow(b1,2)*(pow(LL2,2) - LL2 -1.) + b2);
  als3 = als2 + pow(als0,4)*(pow(b1,3)*(-pow(LL2,3)+5./2.*pow(LL2,2)+2*LL2-1./2.)-3.*b1*b2*LL2 + b3/2.);
  
  return als3;
}


///////////////////////////////////
// evolution from mu=p to mu0=1/a
// Z(mu0)=Z(mu) c(mu0)/c(mu)
// def: c=c(mu)/c(mu0)
// -> Z(mu=1/a) = Z(mu) /c
//////////////////////////////////
double q_evolution_to_RIp_ainv(int Nf,double ainv,double mu_2)
{
  double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
  //mu_2=a2*p2 (adimensional quantity)
  //mu0_2=a2*(1/a^2)=1
  double mu0_2=1;
  
    // alphas @ NNLO
  double alm, al0;
  alm=alphas(Nf,ainv,mu_2)/(4*M_PI);
  al0=alphas(Nf,ainv,mu0_2)/(4*M_PI);
  
  ////////////////////////////////
  // N3LO FORMULA
  // Assuming landau gauge
  ///////////////////////////////////
  if(Nf==2){
    cmu = 1. + 2.03448 * alm + 35.9579 * pow(alm,2) + 1199.16 * pow(alm,3);
    cmu0 = 1. + 2.03448 * al0 + 35.9579 * pow(al0,2) + 1199.16 * pow(al0,3);
  }if(Nf==0){
    cmu = 1. + 2.0303 * alm + 42.1268 * pow(alm,2) + 1728.43 * pow(alm,3);
    cmu0 = 1. + 2.0303 * al0 + 42.1268 * pow(al0,2) + 1728.43 * pow(al0,3);
  }if(Nf==4){
    cmu = 1. + 2.4000 * alm + 29.6724 * pow(alm,2) + 719.141 * pow(alm,3);
    cmu0 = 1. + 2.4000 * al0 + 29.6724 * pow(al0,2) + 719.141 * pow(al0,3);
  }
  
  return cmu/cmu0;
}

double S_evolution_to_RIp_ainv(int Nf,double ainv,double mu_2)
{
  double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
  //mu_2=a2*p2 (adimensional quantity)
  //mu0_2=a2*(1/a^2)=1
  double mu0_2=1;
  
  
  // alphas @ NNLO
  double alm, al0;
  alm=alphas(Nf,ainv,mu_2)/(4*M_PI);
  al0=alphas(Nf,ainv,mu0_2)/(4*M_PI);
  
  
  ////////////////////////////////
  // N3LO FORMULA
  // Assuming landau gauge
  ///////////////////////////////////
  
  
  if(Nf==2){
    cmu = pow(alm,-12./29) * (1. - 8.55727 * alm - 125.423 * pow(alm,2) -
			      3797.71 * pow(alm,3));
    
    
    cmu0 = pow(al0,-12./29) * (1. - 8.55727 * al0 - 125.423 * pow(al0,2) -
			       3797.71 * pow(al0,3));
  }if(Nf==0){
    cmu = pow(alm,-4./11) * (1. - 8.08264 * alm - 151.012 * pow(alm,2) -
			     5247.93 * pow(alm,3));
    
    
    cmu0 = pow(al0,-4./11) * (1. - 8.08264 * al0 - 151.012 * pow(al0,2) -
			      5247.93 * pow(al0,3));
  }if(Nf==4){
    cmu = pow(alm,-12./25) * (1. - 9.38987 * alm - 96.2883 * pow(alm,2) -
			      2403.82 * pow(alm,3));
    
    
    cmu0 = pow(al0,-12./25) * (1. - 9.38987 * al0 - 96.2883 * pow(al0,2) -
			       2403.82 * pow(al0,3));
  }
  
  return cmu/cmu0;
}

double P_evolution_to_RIp_ainv(int Nf,double ainv,double mu_2)
{
    double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
    //mu_2=a2*p2 (adimensional quantity)
    //mu0_2=a2*(1/a^2)=1
    double mu0_2=1;
    
    // alphas @ NNLO
    double alm, al0;
    alm=alphas(Nf,ainv,mu_2)/(4*M_PI);
    al0=alphas(Nf,ainv,mu0_2)/(4*M_PI);
    
    ////////////////////////////////
    // N3LO FORMULA
    // Assuming landau gauge
    ///////////////////////////////////
    if(Nf==2){
      cmu = pow(alm,-12./29) * (1. - 8.55727 * alm - 125.423 * pow(alm,2) -
				3797.71 * pow(alm,3));
      
      cmu0 = pow(al0,-12./29) * (1. - 8.55727 * al0 - 125.423 * pow(al0,2) -
				 3797.71 * pow(al0,3));
    }if(Nf==0){
      cmu = pow(alm,-4./11) * (1. - 8.08264 * alm - 151.012 * pow(alm,2) -
			       5247.93 * pow(alm,3));
      
      
      cmu0 = pow(al0,-4./11) * (1. - 8.08264 * al0 - 151.012 * pow(al0,2) -
				5247.93 * pow(al0,3));
    }if(Nf==4){
      cmu = pow(alm,-12./25) * (1. - 9.38987 * alm - 96.2883 * pow(alm,2) -
				2403.82 * pow(alm,3));
      
      
      cmu0 = pow(al0,-12./25) * (1. - 9.38987 * al0 - 96.2883 * pow(al0,2) -
				 2403.82 * pow(al0,3));
    }
    
    return cmu/cmu0;
}

double T_evolution_to_RIp_ainv(int Nf,double ainv,double mu_2)
{
  double cmu=0.0, cmu0=0.0; // c=cmu/cmu0
  //mu_2=a2*p2 (adimensional quantity)
  //mu0_2=a2*(1/a^2)=1
  double mu0_2=1;
  
  // alphas @ NNLO
  double alm, al0;
  alm=alphas(Nf,ainv,mu_2)/(4*M_PI);
  al0=alphas(Nf,ainv,mu0_2)/(4*M_PI);
  
  ////////////////////////////////
  // N2LO FORMULA
  // Assuming landau gauge
  ///////////////////////////////////
  
  
  if(Nf==2){
    cmu = pow(alm,4./29) * (1. + 2.66852 * alm + 47.9701 * pow(alm,2));
    cmu0 = pow(al0,4./29) * (1. + 2.66852 * al0 + 47.9701 * pow(al0,2));
  }if(Nf==0){
    cmu = pow(alm,4./33) * (1. + 2.53260 * alm + 57.8740 * pow(alm,2));
    cmu0 = pow(al0,4./33) * (1. + 2.53260 * al0 + 57.8740 * pow(al0,2));
  }
  if(Nf==4)
    {
      cmu = pow(alm,4./25) * (1. + 2.91662 * alm + 37.9471 * pow(alm,2));
      cmu0 = pow(al0,4./25) * (1. + 2.91662 * al0 + 37.9471 * pow(al0,2));
  }
  
  return cmu/cmu0;
}
