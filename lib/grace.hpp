#ifndef _GRACE_HPP
#define _GRACE_HPP

#include <fstream>
#include <ave_err.hpp>
#include <tools.hpp>

using namespace std;

namespace grace
{
  enum{WHITE,BLACK,RED,GREEN,BLUE,YELLOW,BROWN,GREY,VIOLET,CYAN,MAGENTA,ORANGE,INDIGO,MAROON,TURQUOISE,GREEN4};
  enum{NO_SYMBOL,CIRCLE,SQUARE,DIAMOND,TRIUP,TRILEFT,TRIDOWN,TRIRIGHT,PLUS,X,STAR};
  enum{NO_LINE,STRAIGHT_LINE};
};

//! class to write a grace file
class grace_file_t : public ofstream
{
  size_t iset;
public:
  //! default constructor
  grace_file_t(const string &path) : ofstream(path),iset(0) {if(!this->good()) CRASH("Unable to open grace file %s",path.c_str());}
  
  //! shift iset
  void shift_iset(size_t how_many=1)
  {iset+=how_many;}
  
  //! start a new set
  void new_set()
  {
    (*this)<<"&"<<endl;
    shift_iset();
  }
  
  //! set a property
  void set_prop(string what){(*this)<<"@s"<<iset<<" "<<what<<endl;}
  
  //! set line style
  void set_line_style(size_t how){set_prop("line type "+to_string(how));}
  void no_line(){set_line_style(grace::NO_LINE);}
  
  //! set colors
  void set_symbol_color(int col){set_prop("symbol color "+to_string(col));}
  void set_symbol_fill_color(int col){set_prop("symbol fill color "+to_string(col));}
  void set_line_color(int col){set_prop("line color "+to_string(col));}
  void set_errorbar_color(int col){set_prop("errorbar color "+to_string(col));}
  void set_colors(int col)
  {
    set_symbol_color(col);
    set_symbol_fill_color(col);
    set_line_color(col);
    set_errorbar_color(col);
  }
  
  //! set symbols
  void set_symbol(int sym){set_prop("symbol "+to_string(sym));}
  void set_no_symbol(){set_symbol(grace::NO_SYMBOL);}
  
  //! form a closed polygon
  void closed_polygon(int fill_col)
  {
    set_prop("fill type 1");
    set_prop("fill color "+to_string(fill_col));
    set_line_color(fill_col);
    set_line_style(grace::STRAIGHT_LINE);
  }
  
  //! write a polygon
  template <class fun_t> void write_polygon(fun_t fun,double xmin,double xmax,size_t npoints,int col=grace::BLACK)
  {
    if(npoints==0) CRASH("NPoints must be different from 0");
    //! x coordinate
    vector<double> x(npoints);
    //! y coordinate
    vec_ave_err_t y(npoints);
    
    //! interval between points
    double dx=(xmax-xmin)/(npoints-1);
    //set x and y
    for(size_t ipoint=0;ipoint<npoints;ipoint++)
      {
	x[ipoint]=xmin+dx*ipoint;
	y[ipoint]=fun(x[ipoint]).ave_err();
      }
    
    //mark a closed polygon
    this->closed_polygon(col);
    //write forward and backward
    for(size_t ipoint=0;ipoint<npoints;ipoint++)         (*this)<<x[ipoint]<<" "<<y[ipoint].ave_minus_err()<<endl;
    for(size_t ipoint=npoints-1;ipoint<npoints;ipoint--) (*this)<<x[ipoint]<<" "<<y[ipoint].ave_plus_err()<<endl;
    
    //create a new set
    this->new_set();
  }
  
  //! write a constant band
  template <class T> void write_constant_band(int xmin,int xmax,const T &c,int col=grace::BLACK)
  {this->write_polygon([&c](double x) -> T {return c;},xmin,xmax,2,col);}
};

//! write a vector of average and error
grace_file_t &operator<<(grace_file_t &out,const vec_ave_err_t &data);

#endif
