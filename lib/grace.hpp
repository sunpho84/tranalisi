#ifndef _GRACE_HPP
#define _GRACE_HPP

#ifndef EXTERN_GRACE
 #define EXTERN_GRACE extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) =__VA_ARGS__
#endif

#include <fstream>
#include <functional>

#include <ave_err.hpp>
#include <functions.hpp>
#include <tools.hpp>

using namespace std;
using namespace std::placeholders;

namespace grace
{
  enum color_t{WHITE,BLACK,RED,GREEN,BLUE,YELLOW,BROWN,GREY,VIOLET,CYAN,MAGENTA,ORANGE,INDIGO,MAROON,TURQUOISE,GREEN4};
  enum symbol_t{NO_SYMBOL,CIRCLE,SQUARE,DIAMOND,TRIUP,TRILEFT,TRIDOWN,TRIRIGHT,PLUS,X,STAR};
  enum line_type_t{NO_LINE_TYPE,STRAIGHT_LINE};
  enum line_style_t{NO_LINE,CONTINUOUS_LINE,SHORT_DASHED_LINE,DASHED_LINE};
  enum fill_type_t{NO_FILL,AS_POLYGON,TO_BASELINE};
  enum symbol_fill_pattern_t{EMPTY_SYMBOL,FILLED_SYMBOL};
  enum settype_t{XY,XYDY,XYDX,XYDXDY};
  
  EXTERN_GRACE symbol_t default_symbol INIT_TO(SQUARE);
  EXTERN_GRACE double default_symbol_size INIT_TO(0.8);
  EXTERN_GRACE double default_symbol_fill_pattern INIT_TO(grace::EMPTY_SYMBOL);
  EXTERN_GRACE color_t default_colour INIT_TO(RED);
  EXTERN_GRACE double default_widths INIT_TO(2);
  EXTERN_GRACE double default_label_size INIT_TO(1.5);
  EXTERN_GRACE vector<grace::color_t> default_color_scheme INIT_TO({grace::RED,grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET});
  EXTERN_GRACE vector<grace::color_t> default_line_color_scheme INIT_TO({grace::RED,grace::BLUE,grace::GREEN4,grace::VIOLET});
  EXTERN_GRACE vector<grace::symbol_t> default_symbol_scheme INIT_TO({grace::TRIUP,grace::CIRCLE,grace::SQUARE,grace::TRIDOWN,grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::STAR});
};

//! class to write a grace file
struct grace_file_t : private ofstream
{
  bool need_close_set;
  
  vector<grace::color_t> color_scheme;
  vector<grace::color_t> line_color_scheme;
  vector<grace::symbol_t> symbol_scheme;
  
  //! get a property
  template <class T>
  T get(const vector<T> &list,const size_t i)
  {
    T out=list[i%list.size()];
    return out;
  }
  
  size_t iset; //!< set id
  string legend; //! legend of the set
  string comment; //! legend of the set
  
  grace::settype_t settype; //!< set type
  grace::symbol_t symbol; //!< symbol
  grace::color_t symbol_color; //!< color for symbol
  bool symbol_fill_pattern; //!< epmty or filled symbols
  grace::color_t symbol_fill_color; //!< color for symbol filling
  double symbol_size; //!< size of the symbol
  double symbol_linewidth; //!< width of the symbol line
  
  grace::line_type_t line_type; //!< no line, straight
  grace::line_style_t line_linestyle; //!< empty or continuous, or dashed in various way
  grace::color_t line_color; //!< color for line
  double linewidth; //!< width of the line
  
  grace::fill_type_t fill_type; //!< filling type
  grace::color_t fill_color; //!< filling type
  
  grace::color_t errorbar_color; //!< color for errorbar
  double errorbar_size; //!< size of the errorbar
  double errorbar_linewidth; //!< width of the bar of therror
  double errorbar_riser_linewidth; //!< width of the bar of the error-riser
  
  string title; //!< title of the plot
  string subtitle; //!< subtitle of the plot
  double title_size; //!< size of the title
  double subtitle_size; //!< size of the subtitle
  bool xaxis_logscale{false}; //!< Whether or not the x-axis is on log scale
  bool yaxis_logscale{false}; //!< Whether or not the y-axis is on log scale
  string xaxis_label; //!< label of the x-axis
  string yaxis_label; //!< lable of the y-axis
  double xaxis_min,xaxis_max; //!< min and max for x-axis
  double yaxis_min,yaxis_max; //!< min and max for y-axis
  double xaxis_label_size; //!< size of the font of the x-axis
  double yaxis_label_size; //!< size of the font of the y-axis
  
  unsigned short int transparency; //!< transparency of all fillings
  
  //! write all props and start a new set
  void close_cur_set();
  
  using ofstream::flush;
  
  size_t cur_col; //!< current color for set (auto-incremented)
  grace::color_t get_col_no_increment();
  
  size_t cur_poly_col; //!< current color for polygon
  grace::color_t get_poly_col_no_increment();
  
  size_t cur_line_col; //!< current color for line
  grace::color_t get_line_col_no_increment();
  
  size_t cur_symbol; //!< current symbol
  grace::symbol_t get_symbol_no_increment();
  
  //! Detect if open
  using ofstream::is_open;
  
  //! direct write
  template <class T>
  friend grace_file_t& operator<<(grace_file_t &file,const T& out)
  {
    static_cast<ofstream&>(file)<<out;
    return file;
  }
  
  //! return col
  size_t get_col() const;
  
  //! return the color scheme
  vector<grace::color_t> get_color_scheme();
  
  //! close current set
  bool get_need_close_set() const;
  
  //! mark the current set as needing close
  void set_need_close_set();
  
  //! set a color scheme
  void set_color_scheme(const initializer_list<grace::color_t> &oth);
  
  void reset_cur_col();
  
  //! set a line color scheme
  void set_line_color_scheme(const initializer_list<grace::color_t> &oth);
  
  void reset_cur_line_col();
  
  //! set a symbol scheme
  void set_symbol_scheme(const initializer_list<grace::symbol_t> &oth);
  
  //! reset all props after creating a new set
  void reset_props();
  
  //! default constructor
  grace_file_t(const string &path="");
  
  //! open
  void open(const string &path);
  
  //! set title of the graph
  void set_title(string label);
  
  //! set subtitle of the graph
  void set_subtitle(string label);
  
  //! set the size of the title of the graph
  void set_title_size(double size);
  
  //! set the size of the subtitle of the graph
  void set_subtitle_size(double size);
  
  //! set logscale of x-axis
  void set_xaxis_logscale(bool yes=true);
  
  //! set title of x-axis
  void set_xaxis_label(string label);
  
  //! set size of the label of the x-axis
  void set_xaxis_label_size(double size);
  
  //! set title of y-axis
  void set_yaxis_label(string label);
  
  //! set size of the label of the y-axis
  void set_yaxis_label_size(double size);
  
  //! set the size of all labels
  void set_all_axis_label_size(double size);
  
  //! set min for x-axis
  void set_xaxis_min(double xmin);
  
  //! set min for y-axis
  void set_yaxis_min(double ymin);
  
  //! set max for x-axis
  void set_xaxis_max(double xmax);
  
  //! set logscale of y-axis
  void set_yaxis_logscale(bool yes=true);
  
  //! set max for y-axis
  void set_yaxis_max(double ymax);
  
  //! set both min and max for x-axis
  void set_xaxis_min_max(double xmin,double xmax);
  
  //! set both min and max for y-axis
  void set_yaxis_min_max(double ymin,double ymax);
  
  //! set the legend of current set
  void set_legend(const string &ext_legend);
  
  //! set the comment of current set
  void set_comment(const string &ext_comment);
  
  //! shift iset
  void shift_iset(size_t how_many=1);
  
  //! start a new set of data type
  void new_data_set(grace::color_t col,grace::symbol_t sym);
  
  //! no color and symbol given
  void new_data_set();
  
  //! set a property
  void write_prop(string what);
  
  //! symbol type
  void set_settype(grace::settype_t ext_settype);
  
  //! set line style and filling
  void set_line_type(grace::line_type_t lt);
  void set_line_style(grace::line_style_t ls);
  void set_no_line();
  void set_fill_type(grace::fill_type_t ft);
  
  //! set colors
  void set_symbol_color(grace::color_t col);
  void set_symbol_fill_color(grace::color_t col);
  void set_symbol_fill_pattern(grace::symbol_fill_pattern_t pat);
  void set_line_color(grace::color_t col);
  void set_fill_color(grace::color_t col);
  void set_errorbar_color(grace::color_t col);
  void set_all_colors(grace::color_t col);
  
  //! set all widths altogether
  void set_all_widths(double w);
  
  //! set the transparency
  void set_transparency(double f);
  
  //! set symbols
  void set_symbol(grace::symbol_t sym);
  void set_no_symbol();
  
  //! form a closed polygon
  void closed_polygon(grace::color_t col);
  
  //! write a polygon
  void write_polygon(const vector<double> &x,const vec_ave_err_t &y,grace::color_t col);
  
  //! write a polygon from a function
  template <class fun_t>
  void write_polygon(const fun_t &fun,double xmin,double xmax,grace::color_t col,size_t npoints=101)
  {
    //! x coordinates
    vector<double> x=vector_grid(xmin,xmax,npoints);
    if(x.size()!=npoints) CRASH("error, x has size %zu and expected %zu points",x.size(),npoints);
    //! y coordinate
    vec_ave_err_t y(npoints);
#pragma omp parallel for
    for(size_t ipoint=0;ipoint<npoints;ipoint++)
      y[ipoint]=fun(x[ipoint]).ave_err();
    
    write_polygon(x,y,col);
  }
  
  //! write a polygon setting automatically col
  void write_polygon(const vector<double> &x,const vec_ave_err_t &y);
  
  //! write a polygon setting automatically col
  template <class fun_t>
  void write_polygon(const fun_t &fun,double xmin,double xmax,size_t npoints=101)
  {
    write_polygon(fun,xmin,xmax,get_poly_col_no_increment(),npoints);
  }
  
  //! mark as a continuos line
  void continuous_line(grace::color_t col);
  
  //! write a line
  template <typename F>
  void write_line(const F &fun,double xmin,double xmax,grace::color_t col,size_t npoints=101)
  {
    close_cur_set();
    
    if(npoints==0) CRASH("NPoints must be different from 0");
    
    //mark a continuous line
    this->continuous_line(col);
    
    double x[npoints],y[npoints];
    
#pragma omp parallel for
    for(size_t ipoint=0;ipoint<npoints;ipoint++)
      {
	if(this->xaxis_logscale)
	  x[ipoint]=xmin*pow(xmax/xmin,ipoint/(double)(npoints-1));
	else
	  x[ipoint]=xmin+(xmax-xmin)/(npoints-1)*ipoint;
	y[ipoint]=fun(x[ipoint]);
      }
    for(size_t ipoint=0;ipoint<npoints;ipoint++)
      (*this)<<x[ipoint]<<" "<<y[ipoint]<<endl;
    need_close_set=true;
  }
  
  template <typename F>
  void write_line(const F &fun,double xmin,double xmax,size_t npoints=101)
  {
    write_line(fun,xmin,xmax,get_line_col_no_increment(),npoints);
  }
  
  //! write a constant band
  template <class T>
  void write_constant_band(double xmin,double xmax,const T &c,grace::color_t col)
  {
    this->write_polygon([&c](double x) -> T {return c;},xmin,xmax,col,2);
  }
  
  template <class T>
  void write_constant_band(double xmin,double xmax,const T &c)
  {
    write_constant_band(xmin,xmax,c,get_line_col_no_increment());
  }
  
  //write a vector of data
  template <class TV,class=enable_if_vector_of_double<TV>>
  void write_vec_ave_err(const TV &x,const vec_ave_err_t &data,grace::color_t col,grace::symbol_t sym)
  {
    new_data_set(col,sym);
    for(size_t i=0;i<data.size();i++)
      if(!std::isnan(data[i].err()))
	(*this)<<x[i]<<" "<<data[i]<<endl;
    set_need_close_set();
  }
  
  void write_vec_ave_err(const vec_ave_err_t &data,grace::color_t col,grace::symbol_t sym);
  
  void write_vec_ave_err(const vec_ave_err_t &data);
  
  template <class TV,class=enable_if_vector_of_double<TV>>
  void write_vec_ave_err(const TV &x,const vec_ave_err_t &data)
  {
    write_vec_ave_err(x,data,get_col_no_increment(),get_symbol_no_increment());
  }
  
  //! write a vector of data
  template <class TV,class=enable_if_vector_of_double<TV>>
  void write_vec_ave_err(const TV &x,const ave_err_t &data,grace::color_t col,grace::symbol_t sym)
  {
    new_data_set(col,sym);
    (*this)<<x<<" "<<data<<endl;
    set_need_close_set();
  }
  
  //! write x and y
  void write_xy(const double x,const double y);
  
  //! write a single data without error on x
  void write_ave_err(const double x,const ave_err_t &data);
  
  //! write a single data with error
  void write_ave_err(const ave_err_t x,const ave_err_t &data);
  
  //! close the file
  void close();
  
  //! destruct
  ~grace_file_t();
};

//! prepare a plot with a band
template <class TV,class T=typename TV::base_type,class fun_t>
void write_fit_plot(grace_file_t &out,double xmin,double xmax,const fun_t &fun,const vector<double> &x,const TV &y)
{
  //(path);
  out.write_polygon(fun,xmin,xmax);
  out.new_data_set();
  out.write_vec_ave_err(x,y.ave_err());
  out.new_data_set();
}
template <class TV,class T=typename TV::base_type,class fun_t>
void write_fit_plot(const string &path,double xmin,double xmax,const fun_t &fun,const vector<double> &x,const TV &y)
{
  grace_file_t out(path);
  write_fit_plot(out,xmin,xmax,fun,x,y);
}

template <class TV,class fun_t,class T=typename TV::base_type>
void write_fit_plot(const string &path,double xmin,double xmax,const fun_t &fun,const TV &y)
{write_fit_plot(path,xmin,xmax,fun,vector_up_to<double>(y.size()),y);}

//! prepare a plot with a polynomial
template <class TV,class T=typename TV::base_type>
void write_poly_fit_plot(const string &path,double xmin,double xmax,const TV &res,const vector<double> &x,const TV &y)
{write_fit_plot(path,xmin,xmax,bind(poly_eval<TV,double>,res,_1),x,y);}
template <class TV,class T=typename TV::base_type> void write_poly_fit_plot(const string &path,double xmin,double xmax,const T&c,const TV &y)
{write_poly_fit_plot(path,xmin,xmax,c,vector_up_to<double>(y.size()),y);}

//! prepare a plot with a polynomial
template <class TV,class T=typename TV::base_type>
void write_constant_fit_plot(const string &path,double xmin,double xmax,const T&c,const vector<double> &x,const TV &y)
{write_fit_plot(path,xmin,xmax,[&c](double x){return c;},x,y);}
template <class TV,class T=typename TV::base_type>
void write_constant_fit_plot(const string &path,double xmin,double xmax,const T&c,const TV &y)
{write_constant_fit_plot(path,xmin,xmax,c,vector_up_to<double>(y.size()),y);}

#undef INIT_TO
#undef EXTERN_GRACE

#endif
