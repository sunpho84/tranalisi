#ifndef _DATA_HPP
#define _DATA_HPP

#include <string>
#include <vector>

#include "node.hpp"

using namespace std;

node_t* parseFormula(const string& formulaString);

/// Data type
struct data_t
{
  const string fileData;
  const double xmin,xmax;
  const node_t* formula;
  const string filePlot;
  
  /// Single data point
  struct point_t
  {
    /// Value of x variable
    double x;
    
    /// Value of y variable
    djack_t y;
    
    /// Error
    double e;
  };
  
  /// Data points
  vector<point_t> points;
  
  /// Constructor
  data_t(const string& fileData,const double xmin,const double xmax,const string& formulaString,const string& filePlot) :
    fileData(fileData),xmin(xmin),xmax(xmax),formula(parseFormula(formulaString)),filePlot(filePlot)
  {
    static int seed=2424;
    
    /// Open the file
    obs_file_t fin(fileData,3,{0,1,2});
    
    /// Count the lines
    const size_t nLines=fin.length(1)/3;
    points.resize(nLines);
    
    for(size_t iLine=0;iLine<nLines;iLine++)
       {
	/// Single line
	const vector<double> l=fin.read(1);
	  points[iLine].x=l[0];
	  points[iLine].y.fill_gauss(l[1],l[2],seed++);
	  points[iLine].e=l[2];
      }
  }
  
  /// Forbids copy constructor
  data_t(const data_t&)=delete;
  
  /// Move constructor
  data_t(data_t&& oth) : fileData(oth.fileData),xmin(oth.xmin),xmax(oth.xmax),formula(oth.formula),filePlot(oth.filePlot)
  {
    oth.formula=nullptr;
    swap(points,oth.points);
  }
  
  /// Destructor
  ~data_t()
  {
    delete formula;
  }
};

#endif
