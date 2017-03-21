#ifndef _GRACE_HPP
#define _GRACE_HPP

#ifndef EXTERN_GRACE
 #define EXTERN_GRACE extern
 #define INIT_TO(...)
#else
 #define INIT_TO(...) =__VA_ARGS__
#endif

#include <fstream>
#include <ave_err.hpp>
#include <tools.hpp>

using namespace std;

namespace grace
{
  enum color_t{WHITE,BLACK,RED,GREEN,BLUE,YELLOW,BROWN,GREY,VIOLET,CYAN,MAGENTA,ORANGE,INDIGO,MAROON,TURQUOISE,GREEN4};
  enum symbol_t{NO_SYMBOL,CIRCLE,SQUARE,DIAMOND,TRIUP,TRILEFT,TRIDOWN,TRIRIGHT,PLUS,X,STAR};
  enum line_type_t{NO_LINE_TYPE,STRAIGHT_LINE};
  enum line_style_t{NO_LINE,CONTINUOUS_LINE,SHORT_DASHED_LINE,DASHED_LINE};
  enum fill_type_t{NO_FILL,AS_POLYGON,TO_BASELINE};
  enum symbol_fill_pattern_t{EMPTY_SYMBOL,FILLED_SYMBOL};
  enum settype_t{XY,XYDY,XYDXDY};
  
  EXTERN_GRACE symbol_t default_symbol INIT_TO(SQUARE);
  EXTERN_GRACE double default_symbol_size INIT_TO(0.8);
  EXTERN_GRACE double default_symbol_fill_pattern INIT_TO(grace::EMPTY_SYMBOL);
  EXTERN_GRACE color_t default_colour INIT_TO(RED);
  EXTERN_GRACE double default_widths INIT_TO(2);
  EXTERN_GRACE double default_label_size INIT_TO(1.5);
  EXTERN_GRACE vector<grace::color_t> default_color_scheme INIT_TO({grace::RED,grace::RED,grace::RED,grace::BLUE,grace::BLUE,grace::GREEN4,grace::VIOLET});
  EXTERN_GRACE vector<grace::color_t> default_line_color_scheme INIT_TO({grace::RED,grace::BLUE,grace::GREEN4});
  EXTERN_GRACE vector<grace::symbol_t> default_symbol_scheme INIT_TO({grace::CIRCLE,grace::SQUARE,grace::DIAMOND,grace::SQUARE,grace::DIAMOND,grace::DIAMOND});
};

//! class to write a grace file
class grace_file_t : public ofstream
{
  bool need_close_set;
  
  vector<grace::color_t> color_scheme;
  vector<grace::color_t> line_color_scheme;
  vector<grace::symbol_t> symbol_scheme;
  //! get a property and increment it
  template <class T> T get_and_increment(const vector<T> &list,size_t &i)
  {
    T out=list[i];
    i=(i+1)%list.size();
    return out;
  }
  
  size_t cur_col; //<! current color for set (auto-incremented)
  grace::color_t get_col_and_increment()
  {return get_and_increment(color_scheme,cur_col);}
  size_t cur_poly_col; //<! current color for polygon
  grace::color_t get_poly_col_and_increment()
  {return get_and_increment(color_scheme,cur_poly_col);}
  size_t cur_line_col; //<! current color for line
  grace::color_t get_line_col_and_increment()
  {return get_and_increment(line_color_scheme,cur_line_col);}
  size_t cur_symbol; //<! current symbol
  grace::symbol_t get_symbol_and_increment()
  {return get_and_increment(symbol_scheme,cur_symbol);}
  
  size_t iset; //<! set id
  string legend; //! legend of the set
  
  grace::symbol_t symbol; //<! symbol
  grace::color_t symbol_color; //<! color for symbol
  bool symbol_fill_pattern; //<! epmty or filled symbols
  grace::color_t symbol_fill_color; //<! color for symbol filling
  double symbol_size; //<! size of the symbol
  double symbol_linewidth; //<! width of the symbol line
  
  grace::line_type_t line_type; //<! no line, straight
  grace::line_style_t line_linestyle; //<! empty or continuous, or dashed in various way
  grace::color_t line_color; //<! color for line
  double linewidth; //<! width of the line
  
  grace::fill_type_t fill_type; //<! filling type
  grace::color_t fill_color; //<! filling type
  
  grace::color_t errorbar_color; //<! color for errorbar
  double errorbar_size; //<! size of the errorbar
  double errorbar_linewidth; //<! width of the bar of therror
  double errorbar_riser_linewidth; //<! width of the bar of the error-riser
  
  string title; //<! title of the plot
  string subtitle; //<! subtitle of the plot
  double title_size; //! size of the title
  double subtitle_size; //! size of the subtitle
  string xaxis_label; //<! label of the x-axis
  string yaxis_label; //<! lable of the y-axis
  double xaxis_min,xaxis_max; //<! min and max for x-axis
  double yaxis_min,yaxis_max; //<! min and max for y-axis
  double xaxis_label_size; //<! size of the font of the x-axis
  double yaxis_label_size; //<! size of the font of the y-axis
  
  unsigned short int transparency; //<! transparency of all fillings
  
public:
  //! set a color scheme
  void set_color_scheme(const initializer_list<grace::color_t> &oth)
  {color_scheme.assign(oth.begin(),oth.end());}
  void reset_cur_col()
  {cur_col=0;}

  //! set a line color scheme
  void set_line_color_scheme(const initializer_list<grace::color_t> &oth)
  {line_color_scheme.assign(oth.begin(),oth.end());}
  void reset_cur_line_col()
  {cur_line_col=0;}
  
  //! set a symbol scheme
  void set_symbol_scheme(const initializer_list<grace::symbol_t> &oth)
  {symbol_scheme.assign(oth.begin(),oth.end());}
  
  //! reset all props after creating a new set
  void reset_props()
  {
    line_type=grace::STRAIGHT_LINE;
    line_linestyle=grace::CONTINUOUS_LINE;
    fill_type=grace::NO_FILL;
    symbol=grace::default_symbol;
    symbol_size=grace::default_symbol_size;
    symbol_fill_pattern=grace::default_symbol_fill_pattern;
    set_all_colors(grace::default_colour);
    transparency=255;
    errorbar_size=0.5;
    set_all_widths(grace::default_widths);
  }
  
  //! default constructor
  grace_file_t(const string &path) :
    ofstream(path),
    need_close_set(false),
    color_scheme(grace::default_color_scheme),
    line_color_scheme(grace::default_line_color_scheme),
    symbol_scheme(grace::default_symbol_scheme),
    cur_col(0),
    cur_poly_col(0),
    cur_line_col(0),
    cur_symbol(0),
    iset(0),
    xaxis_min(0),
    xaxis_max(1),
    yaxis_min(0),
    yaxis_max(1)
  {
    title_size=
      subtitle_size=
      xaxis_label_size=
      yaxis_label_size=
      grace::default_label_size;
    
    reset_props();
    if(!this->good()) CRASH("Unable to open grace file %s",path.c_str());
  }
  
  //! set title of the graph
  void set_title(string label)
  {title=label;}
  
  //! set subtitle of the graph
  void set_subtitle(string label)
  {subtitle=label;}
  
  //! set the size of the title of the graph
  void set_title_size(double size)
  {title_size=size;}
  
  //! set the size of the subtitle of the graph
  void set_subtitle_size(double size)
  {subtitle_size=size;}
  
  //! set title of x-axis
  void set_xaxis_label(string label)
  {xaxis_label=label;}
  
  //! set size of the label of the x-axis
  void set_xaxis_label_size(double size)
  {xaxis_label_size=size;}
  
  //! set title of y-axis
  void set_yaxis_label(string label)
  {yaxis_label=label;}
  
  //! set size of the label of the y-axis
  void set_yaxis_label_size(double size)
  {yaxis_label_size=size;}
  
  //! set the size of all labels
  void set_all_axis_label_size(double size)
  {
    set_xaxis_label_size(size);
    set_yaxis_label_size(size);
  }
  
  //! set min for x-axis
  void set_xaxis_min(double xmin)
  {xaxis_min=xmin;}
  
  //! set min for y-axis
  void set_yaxis_min(double ymin)
  {yaxis_min=ymin;}
  
  //! set max for x-axis
  void set_xaxis_max(double xmax)
  {xaxis_max=xmax;}
  
  //! set max for y-axis
  void set_yaxis_max(double ymax)
  {yaxis_max=ymax;}
  
  //! set both min and max for x-axis
  void set_xaxis_min_max(double xmin,double xmax)
  {
    set_xaxis_min(xmin);
    set_xaxis_max(xmax);
  }
  
  //! set both min and max for y-axis
  void set_yaxis_min_max(double ymin,double ymax)
  {
    set_yaxis_min(ymin);
    set_yaxis_max(ymax);
  }
  
  //! set the legend of current set
  void set_legend(const string &ext_legend)
  {legend=ext_legend;}
  
  //! shift iset
  void shift_iset(size_t how_many=1)
  {iset+=how_many;}
  
  //! write all props and start a new set
  void close_cur_set()
  {
    if(need_close_set)
      {
	//line props
	write_prop("line type "+to_string(line_type));
	write_prop("linewidth "+to_string(linewidth));
	write_prop("line linestyle "+to_string(line_linestyle));
	write_prop("line color "+to_string(line_color));
	//fill props
	write_prop("fill color "+to_string(fill_color));
	write_prop("fill type "+to_string(fill_type));
	//symbol props
	write_prop("symbol "+to_string(symbol));
	write_prop("symbol size "+to_string(symbol_size));
	write_prop("symbol linewidth "+to_string(symbol_linewidth));
	write_prop("symbol color "+to_string(symbol_color));
	write_prop("symbol fill color "+to_string(symbol_fill_color));
	write_prop("symbol fill pattern "+to_string(symbol_fill_pattern));
	//error props
	write_prop("errorbar color "+to_string(errorbar_color));
	write_prop("errorbar size "+to_string(errorbar_size));
	write_prop("errorbar linewidth "+to_string(errorbar_linewidth));
	write_prop("errorbar riser linewidth "+to_string(errorbar_riser_linewidth));
	//transparency
	(*this)<<"#QTGRACE_ADDITIONAL_PARAMETER: G 0 S "<<iset<<" ALPHA_CHANNELS {"<<transparency<<";"<<transparency<<";255;255;255;255}"<<endl;
	//write legend
	write_prop("legend \""+legend+"\"");
	//point to the next set
	(*this)<<"&"<<endl;
	shift_iset();
	(*this)<<"@target G0.S"<<iset<<endl;
	//reset all props and force next call to ignore
	reset_props();
	need_close_set=false;
      }
  }
  
  //! start a new set of data type
  void new_data_set(grace::color_t col,grace::symbol_t sym)
  {
    close_cur_set();
    
    set_settype(grace::XYDY);
    set_line_style(grace::NO_LINE);
    set_all_colors(col);
    set_symbol(sym);
    
    need_close_set=true;
  }
  void new_data_set() //<! no color and symbol given
  {new_data_set(get_col_and_increment(),get_symbol_and_increment());}
  
  //! set a property
  void write_prop(string what) {(*this)<<"@s"<<iset<<" "<<what<<endl;}
  
  //! symbol type
  void set_settype(grace::settype_t settype)
  {
    string how_s;
    switch(settype)
      {
      case grace::XY: how_s="xy";break;
      case grace::XYDY: how_s="xydy";break;
      case grace::XYDXDY: how_s="xydxdy";break;
      default: CRASH("Unknown type %d",settype);
      }
    (*this)<<"@type "<<how_s<<endl;
  }
  
  //! set line style and filling
  void set_line_type(grace::line_type_t lt) {line_type=lt;}
  void set_line_style(grace::line_style_t ls) {line_linestyle=ls;}
  void set_no_line() {set_line_style(grace::NO_LINE);}
  void set_fill_type(grace::fill_type_t ft) {fill_type=ft;}
  
  //! set colors
  void set_symbol_color(grace::color_t col) {symbol_color=col;}
  void set_symbol_fill_color(grace::color_t col) {symbol_fill_color=col;}
  void set_symbol_fill_pattern(grace::symbol_fill_pattern_t pat) {symbol_fill_pattern=pat;}
  void set_line_color(grace::color_t col) {line_color=col;}
  void set_fill_color(grace::color_t col) {fill_color=col;}
  void set_errorbar_color(grace::color_t col) {errorbar_color=col;}
  void set_all_colors(grace::color_t col)
  {
    set_symbol_color(col);
    set_symbol_fill_color(col);
    set_fill_color(col);
    set_line_color(col);
    set_errorbar_color(col);
  }
  
  //! set all widths altogether
  void set_all_widths(double w)
  {
    linewidth=
      symbol_linewidth=
      errorbar_linewidth=
      errorbar_riser_linewidth=
      w;
  }
  
  //! set the transparency
  void set_transparency(double f)
  {
    if(f>1 or f<0) CRASH("f=%lg, must be within [0;1]",f);
    transparency=(unsigned short int)(f*255);
  }
  
  //! set symbols
  void set_symbol(grace::symbol_t sym) {symbol=sym;}
  void set_no_symbol() {set_symbol(grace::NO_SYMBOL);}
  
  //! form a closed polygon
  void closed_polygon(grace::color_t col)
  {
    set_settype(grace::XY);
    set_fill_type(grace::AS_POLYGON);
    set_no_line();
    set_all_colors(col);
    set_transparency(0.4);
    set_no_symbol();
  }
  
  //! write a polygon
  template <class fun_t> void write_polygon(fun_t &fun,double xmin,double xmax,grace::color_t col,size_t npoints=100)
  {
    close_cur_set();
    
    //mark a closed polygon
    this->closed_polygon(col);
    
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
    
    //write forward and backward
    for(size_t ipoint=0;ipoint<npoints;ipoint++)         (*this)<<x[ipoint]<<" "<<y[ipoint].ave_minus_err()<<endl;
    for(size_t ipoint=npoints-1;ipoint<npoints;ipoint--) (*this)<<x[ipoint]<<" "<<y[ipoint].ave_plus_err()<<endl;
    
    need_close_set=true;
  }
  template <class fun_t> void write_polygon(fun_t &fun,double xmin,double xmax,size_t npoints=100)
  {write_polygon(fun,xmin,xmax,get_poly_col_and_increment(),npoints);}
  
  //! mark as a continuos line
  void continuous_line(grace::color_t col)
  {
    set_settype(grace::XY);
    set_all_colors(col);
    set_no_symbol();
  }
  
  //! write a line
  template <class fun_t> void write_line(fun_t fun,double xmin,double xmax,grace::color_t col,size_t npoints=100)
  {
    close_cur_set();
    
    if(npoints==0) CRASH("NPoints must be different from 0");
    
    //mark a continuous line
    this->continuous_line(col);
    
    for(size_t ipoint=0;ipoint<npoints;ipoint++)
      {
	double x=xmin+(xmax-xmin)/(npoints-1)*ipoint;
	double y=fun(x);
	(*this)<<x<<" "<<y<<endl;
      }
    need_close_set=true;
  }
  template <class fun_t> void write_line(fun_t fun,double xmin,double xmax,size_t npoints=100)
  {write_line(fun,xmin,xmax,get_line_col_and_increment(),npoints);}
  
  //! write a constant band
  template <class T> void write_constant_band(int xmin,int xmax,const T &c,grace::color_t col)
  {this->write_polygon([&c](double x) -> T {return c;},xmin,xmax,col,2);}
  template <class T> void write_constant_band(int xmin,int xmax,const T &c)
  {write_constant_band(xmin,xmax,c,get_col_and_increment());}
  
  //write a vector of data
  void write_vec_ave_err(const vec_ave_err_t &data,grace::color_t col,grace::symbol_t sym)
  {
    new_data_set(col,sym);
    for(size_t i=0;i<data.size();i++)
      if(!std::isnan(data[i].err))
	(*this)<<i<<" "<<data[i]<<endl;
  }
  void write_vec_ave_err(const vec_ave_err_t &data)
  {write_vec_ave_err(data,get_col_and_increment(),get_symbol_and_increment());}
  
  //write a single data
  void write_ave_err(const double x,const ave_err_t &data,grace::color_t col,grace::symbol_t sym)
  {
    new_data_set(col,sym);
    (*this)<<x<<" "<<data<<endl;
  }
  void write_ave_err(const double x,const ave_err_t &data)
  {write_ave_err(x,data,get_col_and_increment(),get_symbol_and_increment());}
  
  //write a single data
  void write_ave_err(const ave_err_t x,const ave_err_t &data,grace::color_t col,grace::symbol_t sym)
  {
    new_data_set(col,sym);
    set_settype(grace::XYDXDY);
    (*this)<<x.ave<<" "<<data.ave<<" "<<x.err<<" "<<data.err<<endl;
  }
  void write_ave_err(const ave_err_t x,const ave_err_t &data)
  {write_ave_err(x,data,get_col_and_increment(),get_symbol_and_increment());}
  
  //! close the file
  void close()
  {
    close_cur_set();
    
    (*this)<<"@title \""<<title<<"\""<<endl;
    (*this)<<"@subtitle \""<<subtitle<<"\""<<endl;
    (*this)<<"@title size "<<title_size<<endl;
    (*this)<<"@subtitle size "<<subtitle_size<<endl;
    (*this)<<"@xaxis label \""<<xaxis_label<<"\""<<endl;
    (*this)<<"@xaxis label char size "<<xaxis_label_size<<endl;
    (*this)<<"@yaxis label \""<<yaxis_label<<"\""<<endl;
    (*this)<<"@yaxis label char size "<<yaxis_label_size<<endl;
    (*this)<<"@world "<<xaxis_min<<","<<yaxis_min<<","<<xaxis_max<<","<<yaxis_max<<endl;
    
    this->ofstream::close();
  }
  
  //! destruct
  ~grace_file_t()
  {close();}
};

//! prepare a plot with a band
template <class TV,class T=typename TV::base_type,class fun_t> void write_fit_plot(const string &path,double xmin,double xmax,fun_t fun,const TV &v)
{
  grace_file_t out(path);
  out.write_polygon(fun,xmin,xmax);
  out.new_data_set();
  out.write_vec_ave_err(v.ave_err());
  out.new_data_set();
}
//! prepare a plot with a constant
template <class TV,class T=typename TV::base_type> void write_constant_fit_plot(const string &path,double xmin,double xmax,const T&c,const TV &v)
{write_fit_plot(path,xmin,xmax,[&c](double x) -> T {return c;},v);}

#undef INIT_TO
#undef EXTERN_GRACE

#endif
