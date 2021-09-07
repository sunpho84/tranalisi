#ifndef _DISTANCEMULTIPLICITY_HPP
#define _DISTANCEMULTIPLICITY_HPP

#include <vector>

using namespace std;

/// Computes and keep track of the multiplicity of points, up da a maximal distance
struct DistanceMultiplicity
{
  //// Maximal computed distance
  const int maxDistance;
  
  /// Square of the maximal distance
  const int maxDistance2;
  
  /// Table of multiplicity
  vector<int> multiplicityTable;
  
  /// Fills the multiplicity table
  void fillMultiplicityTable()
  {
    /// Multiplicity of the point, according to the number of zero coordinates
    constexpr int zeroCoordMultiplicity[4]={8,4,2,1};
    
    /// Multiplicity of the point, according to the number of equal coordinates
    constexpr int equalCoordMultiplicity[4]={6,3,0,1};
    
    multiplicityTable.resize(maxDistance2,0);
    
    for(int i=0;i<=maxDistance;i++)
      for(int j=i;j<=maxDistance;j++)
	for(int k=j;k<=maxDistance;k++)
	  {
	    /// Number of zero coordinates
	    const int nZeroCoords=
	      (i==0)+(j==0)+(k==0);
	    
	    /// Number of euqal coordinates
	    const int nEqualCoords=
	      (i==j)+(i==k)+(j==k);
	    
	    /// Squared distance
	    const int dist2=
	      i*i+j*j+k*k;
	    
	    if(dist2<maxDistance2)
	      multiplicityTable[dist2]+=
		zeroCoordMultiplicity[nZeroCoords]*equalCoordMultiplicity[nEqualCoords];
	  }
  }
  
  /// Constructor
  DistanceMultiplicity(const int& maxDistance=30) :
    maxDistance(maxDistance),
    maxDistance2(maxDistance*maxDistance+1)
  {
    fillMultiplicityTable();
  }
  
  /// Access the table
  const int& operator()(const int& i) const
  {
    return multiplicityTable[i];
  }
};

#endif
