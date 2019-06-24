#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_GRACE
 #include <grace.hpp>

grace::color_t grace_file_t::get_col_no_increment()
{
  return get(color_scheme,cur_col);
}

grace::color_t grace_file_t::get_poly_col_no_increment()
{
  return get(color_scheme,cur_poly_col);
}

grace::color_t grace_file_t::get_line_col_no_increment()
{
  return get(line_color_scheme,cur_line_col);
}

grace::symbol_t grace_file_t::get_symbol_no_increment()
{
  return get(symbol_scheme,cur_symbol);
}

void grace_file_t::close_cur_set()
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
      
      cur_col++;
      cur_poly_col++;
      cur_line_col++;
      cur_symbol++;
    }
}

size_t grace_file_t::get_col() const
{
  return cur_col;
}

vector<grace::color_t> grace_file_t::get_color_scheme()
{
  return color_scheme;
}

bool grace_file_t::get_need_close_set() const
{
  return need_close_set;
}

void grace_file_t::set_need_close_set()
{
  need_close_set=true;
}

void grace_file_t::set_color_scheme(const initializer_list<grace::color_t> &oth)
{
  color_scheme.assign(oth.begin(),oth.end());
}

void grace_file_t::reset_cur_col()
{
  cur_col=0;
}

void grace_file_t::set_line_color_scheme(const initializer_list<grace::color_t> &oth)
{
  line_color_scheme.assign(oth.begin(),oth.end());
}

void grace_file_t::reset_cur_line_col()
{
  cur_line_col=0;
}

void grace_file_t::set_symbol_scheme(const initializer_list<grace::symbol_t> &oth)
{
  symbol_scheme.assign(oth.begin(),oth.end());
}

void grace_file_t::reset_props()
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

grace_file_t::grace_file_t(const string &path) :
  need_close_set(false),
  color_scheme(grace::default_color_scheme),
  line_color_scheme(grace::default_line_color_scheme),
  symbol_scheme(grace::default_symbol_scheme),
  iset(0),
  settype(grace::XY),
  xaxis_min(0),
  xaxis_max(1),
  yaxis_min(0),
  yaxis_max(1),
  cur_col(0),
  cur_poly_col(0),
  cur_line_col(0),
  cur_symbol(0)
{
  if(path!="") this->open(path);
  
  title_size=
    subtitle_size=
    xaxis_label_size=
    yaxis_label_size=
    grace::default_label_size;
  
  this->precision(16);
  
  reset_props();
}

void grace_file_t::open(const string &path)
{
  ofstream::open(path);
  if(!this->good()) CRASH("Unable to open grace file %s",path.c_str());
}

void grace_file_t::set_title(string label)
{
  title=label;
}

void grace_file_t::set_subtitle(string label)
{
  subtitle=label;
}

void grace_file_t::set_title_size(double size)
{
  title_size=size;
}

void grace_file_t::set_subtitle_size(double size)
{
  subtitle_size=size;
}

void grace_file_t::set_xaxis_logscale(bool yes)
{
  xaxis_logscale=yes;
  
  if(yes)
    (*this)<<"@    xaxes scale Logarithmic\n";
  else
    (*this)<<"@    xaxes scale Normal\n";
}

void grace_file_t::set_xaxis_label(string label)
{
  xaxis_label=label;
}

void grace_file_t::set_xaxis_label_size(double size)
{
  xaxis_label_size=size;
}

void grace_file_t::set_yaxis_label(string label)
{
  yaxis_label=label;
}

void grace_file_t::set_yaxis_label_size(double size)
{
  yaxis_label_size=size;
}

void grace_file_t::set_all_axis_label_size(double size)
{
  set_xaxis_label_size(size);
  set_yaxis_label_size(size);
}
  
void grace_file_t::set_xaxis_min(double xmin)
{
  xaxis_min=xmin;
}

void grace_file_t::set_yaxis_min(double ymin)
{
  yaxis_min=ymin;
}

void grace_file_t::set_xaxis_max(double xmax)
{
  xaxis_max=xmax;
}

void grace_file_t::set_yaxis_logscale(bool yes)
{
  yaxis_logscale=yes;
  
  if(yes)
    (*this)<<"@    yaxes scale Logarithmic\n";
  else
    (*this)<<"@    yaxes scale Normal\n";
}

void grace_file_t::set_yaxis_max(double ymax)
{
  yaxis_max=ymax;
}

void grace_file_t::set_xaxis_min_max(double xmin,double xmax)
{
  set_xaxis_min(xmin);
  set_xaxis_max(xmax);
}

void grace_file_t::set_yaxis_min_max(double ymin,double ymax)
{
  set_yaxis_min(ymin);
  set_yaxis_max(ymax);
}

void grace_file_t::set_legend(const string &ext_legend)
{
  legend=ext_legend;
}

void grace_file_t::shift_iset(size_t how_many)
{
  iset+=how_many;
}

void grace_file_t::new_data_set(grace::color_t col,grace::symbol_t sym)
{
  close_cur_set();
  
  set_legend("");
  set_settype(grace::XYDY);
  set_all_colors(col);
  set_symbol(sym);
  
  need_close_set=true;
}

void grace_file_t::new_data_set()
{
  new_data_set(get_col_no_increment(),get_symbol_no_increment());
}

void grace_file_t::write_prop(string what)
{
  (*this)<<"@s"<<iset<<" "<<what<<endl;
}

void grace_file_t::set_settype(grace::settype_t ext_settype)
{
  if(this->settype!=ext_settype)
    {
      string how_s;
      switch(ext_settype)
	{
	case grace::XY: how_s="xy";break;
	case grace::XYDX: how_s="xydx";break;
	case grace::XYDY: how_s="xydy";break;
	case grace::XYDXDY: how_s="xydxdy";break;
	default: CRASH("Unknown type %d",settype);
	}
      (*this)<<"@type "<<how_s<<endl;
    }
  
  this->settype=ext_settype;
}

void grace_file_t::set_line_type(grace::line_type_t lt)
{
  line_type=lt;
}

void grace_file_t::set_line_style(grace::line_style_t ls)
{
  line_linestyle=ls;
}

void grace_file_t::set_no_line()
{
  set_line_style(grace::NO_LINE);
}

void grace_file_t::set_fill_type(grace::fill_type_t ft)
{
  fill_type=ft;
}

void grace_file_t::set_symbol_color(grace::color_t col)
{
  symbol_color=col;
}

void grace_file_t::set_symbol_fill_color(grace::color_t col)
{
  symbol_fill_color=col;
}

void grace_file_t::set_symbol_fill_pattern(grace::symbol_fill_pattern_t pat)
{
  symbol_fill_pattern=pat;
}

void grace_file_t::set_line_color(grace::color_t col)
{
  line_color=col;
}

void grace_file_t::set_fill_color(grace::color_t col)
{
  fill_color=col;
}

void grace_file_t::set_errorbar_color(grace::color_t col)
{
  errorbar_color=col;
}

void grace_file_t::set_all_colors(grace::color_t col)
{
  set_symbol_color(col);
  set_symbol_fill_color(col);
  set_fill_color(col);
  set_line_color(col);
  set_errorbar_color(col);
}

void grace_file_t::set_all_widths(double w)
{
  linewidth=
    symbol_linewidth=
    errorbar_linewidth=
    errorbar_riser_linewidth=
    w;
}

void grace_file_t::set_transparency(double f)
{
  if(f>1 or f<0) CRASH("f=%lg, must be within [0;1]",f);
  transparency=(unsigned short int)(f*255);
}

void grace_file_t::set_symbol(grace::symbol_t sym)
{
  symbol=sym;
}

void grace_file_t::set_no_symbol()
{
  set_symbol(grace::NO_SYMBOL);
}

void grace_file_t::closed_polygon(grace::color_t col)
{
  set_settype(grace::XY);
  set_fill_type(grace::AS_POLYGON);
  set_no_line();
  set_all_colors(col);
  set_transparency(0.4);
  set_no_symbol();
}

void grace_file_t::write_polygon(const vector<double> &x,const vec_ave_err_t &y,grace::color_t col)
{
  size_t npoints=x.size();
  if(y.size()!=npoints) CRASH("x and y have different size, %zu and %zu",npoints,y.size());
  
  //mark a closed polygon
  close_cur_set();
  this->closed_polygon(col);
  
  //write forward and backward
  for(size_t ipoint=0;ipoint<npoints;ipoint++)         (*this)<<x[ipoint]<<" "<<y[ipoint].ave_minus_err()<<endl;
  for(size_t ipoint=npoints-1;ipoint<npoints;ipoint--) (*this)<<x[ipoint]<<" "<<y[ipoint].ave_plus_err()<<endl;
  
  need_close_set=true;
  close_cur_set();
}

void grace_file_t::write_polygon(const vector<double> &x,const vec_ave_err_t &y)
{
  write_polygon(x,y,get_poly_col_no_increment());
}

void grace_file_t::continuous_line(grace::color_t col)
{
  set_settype(grace::XY);
  set_all_colors(col);
  set_no_symbol();
}

void grace_file_t::write_vec_ave_err(const vec_ave_err_t &data,grace::color_t col,grace::symbol_t sym)
{
  write_vec_ave_err(vector_up_to<double>(data.size()),data,col,sym);
}

void grace_file_t::write_vec_ave_err(const vec_ave_err_t &data)
{
  write_vec_ave_err(data,get_col_no_increment(),get_symbol_no_increment());
}

void grace_file_t::write_xy(const double x,const double y)
{
  set_need_close_set();
  set_settype(grace::XY);
  (*this)<<x<<" "<<y<<endl;
}

void grace_file_t::write_ave_err(const double x,const ave_err_t &data)
{
  set_need_close_set();
  set_settype(grace::XYDY);
  (*this)<<x<<" "<<data.ave()<<" "<<data.err()<<endl;
}

void grace_file_t::write_ave_err(const ave_err_t x,const ave_err_t &data)
{
  set_need_close_set();
  set_settype(grace::XYDXDY);
  (*this)<<x.ave()<<" "<<data.ave()<<" "<<x.err()<<" "<<data.err()<<endl;
}

void grace_file_t::close()
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

grace_file_t::~grace_file_t()
{
  close();
}
