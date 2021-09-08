#ifndef _PHIFUNCTION_HPP
#define _PHIFUNCTION_HPP

#include <fit.hpp>

#include <luscherZeta.hpp>

/// Computes the r.h.s of eq.A.3 of Nuclear Physics B364 (1991) 237-251
struct TanPhiCalculator
{
  /// Computes the Zeta function
  const LuscherZetaCalculator luscherZetaCalculator;
  
  /// Approximation according to eq. A.4
  double approximation(const double& q2) const
  {
    /// Equation A.5
    const double x=
      q2*(8.91363-
	  q2*(16.5323-
	      q2*(8.402-
		  q2*6.95)));
    
    /// Numerator
    const double num=
      2*M_PI*M_PI*pow(q2,1.5);
    
    /// Denominator
    const double den=
      1+x;
    
    return
      num/den;
  }
  
  /// Full evalation
  double full(const double& q2) const
  {
    /// Square root of q2
    const double q=
      sqrt(q2);
    
    /// Pi^(3/2)
    static const double pi32=
      pow(M_PI,1.5);
    
    /// Numerator
    const double num=
      pi32*q;
    
    /// Denominator
    const double den=
      luscherZetaCalculator(q2);
    
    return
      -num/den;
  }
  
  /// Calculate tan(phi)
  double operator()(const double& q2) const
  {
    if(q2<=0.01)
      return approximation(q2);
    else
      return full(q2);
  }
};

/// Structure to compute phi
struct PhiCalculator
{
  /// Calculator of the tan
  TanPhiCalculator tanPhiCalculator;
  
  /// Inferior limit
  const double xMin;
  
  /// Superior limit
  const double xMax;
  
  /// Number of points
  const int nPoints;
  
  /// Interval
  const double dX;
  
  /// Numer of known ones (without shifting)
  const vector<double> knownPhiOnes;
  
  /// Temporarily store the points
  vector<double> y;
  
  /// I-th point's abscissa
  double x(const int& iPoint) const
  {
    return
      xMin+dX*iPoint;
  }
  
  /// Finds all points where the function is approximately one in modulo
  vector<double> findPhiOnes() const
  {
    /// List of ones position to be returned
    vector<double> res;
    
    /// Previous value
    double prev=0.0;
    
    for(double x=0.1;x<xMax+2.0;x+=dX)
      {
	///Current value
	const double cur=
	  atan(tanPhiCalculator(x));
	
	if((prev<1.0 and cur>1.0) or (prev<-1.0 and cur>-1.0))
	  {
	    /// Position
	    const double pos=
	      x-dX/2.0;
	    
	    res.push_back(pos);
	  }
	
	prev=cur;
      }
    
    return res;
  }
  
  /// Constructor
  PhiCalculator(const double& xMin,
		const double& xMax,
		const int& nIntervals) :
    xMin(xMin),
    xMax(xMax),
    nPoints(nIntervals+1),
    dX((xMax-xMin)/nIntervals),
    knownPhiOnes(findPhiOnes()),
    y(nPoints)
  {
    cout<<"Ones:"<<endl;
    cout<<knownPhiOnes<<endl;
    
  }
  
  /// Computes putting the correct shift to provide continuity
  double operator()(const double& x) const
  {
    if(x>=xMax)
      CRASH("Cannot compute beyond xMax");
    
    const double phi=
      atan(tanPhiCalculator(x));
    
    int iOnes=0;
    while(iOnes<(int)knownPhiOnes.size()-1 and knownPhiOnes[iOnes+1]<x)
      iOnes++;
    
    int b=
      (iOnes+1)/2+(iOnes%2==0 and phi<0);
      
    const double shiftedPhi=
      phi+M_PI*b;
    
    return
      shiftedPhi;
  }
};

  // /// Previous value of phi
    // double oldPhi=
    //   0.0;
    
    // /// Correction needed to smooth phi
    // int iCorr=0;
    
    // for(int iPoint=0;iPoint<nPoints;iPoint++)
    //   {
    // 	/// Abscissa
    // 	const double x=
    // 	  dX*iPoint+xMin;
	
    // 	/// Tan of phi
    // 	const double tanPhi=
    // 	  tanPhiCalculator(x);
	
    // 	// const double phi=atan(tanPhi);
    // 	const double phi=
    // 	  atan(tanPhi);
	
    // 	if(iPoint and x>0.1 and phi*oldPhi<-1.0)
    // 	  iCorr++;
	
    // 	/// Smoothing correction
    // 	const double corr=
    // 	  phi+iCorr*M_PI;
	
    // 	y[iPoint]=
    // 	  corr/(M_PI*x+1e-16);
	
    // 	oldPhi=phi;
      // }
//   }
  
//   vector<double> getSpline(const int& iX,
// 			   const int& degree) const
//   {
//     vector<double> Al(2*degree+1,0.0);
//     vector<double> c(degree+1,0.0);
  
//     for(int dI=-(degree-1)/2;dI<=(degree+1)/2;dI++)
//       {
// 	const int iPoint=
// 	  iX+dI;
	
// 	const double x=
// 	  dI*dX;
	
//         /// Weight
//         double w=
// 	  1.0;
	
//         for(int f=0;f<=2*degree;f++)
//           {
//             Al[f]+=w;
//             if(f<=degree)
// 	      c[f]+=y[iPoint]*w;
//             w*=x;
//           }
//       }
  
//   vector<double> A((degree+1)*(degree+1));
//   for(int i=0;i<=degree;i++)
//     for(int j=0;j<=degree;j++)
//       A[i*(degree+1)+j]=Al[i+j];
  
//   //
  
//   for(int i=0;i<degree+1;i++)
//     {
//       double C=A[i*(degree+1)+i];
//       for(int j=i;j<degree+1;j++)
// 	A[i*(degree+1)+j]/=C;
//       c[i]/=C;
      
//       for(int k=i+1;k<degree+1;k++)
//         {
//           double C=A[k*(degree+1)+i];
//           for(int j=i;j<degree+1;j++)
// 	    A[k*(degree+1)+j]-=A[i*(degree+1)+j]*C;
//           c[k]-=C*c[i];
//         }
//     }
  
//   vector<double> res(degree+1);
//   for(int k=degree;k>=0;k--)
//     {
//       double S=
// 	0.0;
      
//       for(int i=k+1;i<degree+1;i++)
// 	S+=A[k*(degree+1)+i]*res[i];
//       res[k]=c[k]-S;
//     }
  
//   return
//     res;
//   }
  
//     const int iX=
//       floor((x-xMin)/dX);
    
//     constexpr int degree=
//       3;
    
//     static_assert(degree%2,"degree must be odd");
    
//     const int offset=
//       (degree-1)/2;
    
//     const int i0=
//       std::max(offset,std::min(iX,nPoints-offset-2));
    
//     // poly_fit(const vector<double> &, const TV &y, int d)
    
//     // const int i1=
//     //   i0+1;
    
//     // const double& y0=
//     //   y[i0];
    
//     // return y0;
    
//     // const double& y1=
//     //   y[i1];
    
//     // const double& x1=
//     //   i1*dX+xMin;
    
//     // const double a=
//     //   (y1-y0)/dX;
    
//     // const double b=
//     //   y1-x1*a;
    
//     const vector<double> coeffs=
//       getSpline(i0,degree);
    
//     double res=
//       0.0;
    
//     double rx=
//       1.0;
    
//     const double dist=
//       x-(i0*dX+xMin);
    
//     for(int i=0;i<=degree;i++)
//       {
// 	res+=coeffs[i]*rx;
// 	rx*=dist;
//       }
    
//     return
//       res;
//   }
// };

#endif
