#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <ave_err.hpp>
#include <grace.hpp>
#include <oper.hpp>

//! write to a path
void vec_ave_err_t::write(const vector<double> &x,const string &path) const
{
  grace_file_t gr(path);
  gr.write_vec_ave_err(x,*this);
}

//! shift by 10 left(if positive) or right(if negative)
double shift_by_10(double x,int digit)
{
  double fact=(digit>0)?10.0:0.1;
  for(int i=0;i<abs(digit);i++) x*=fact;
  return x;
}

//! return x rounded to the digit: (18.32,-1)=18.3, (18.32,1)=20
double round_to_digit(double x,int digit)
{return shift_by_10(floor(shift_by_10(x,-digit)+0.5),digit);}

//! print the average and errors considering all rounding
string smart_print(double ave,const vector<double> &errors,int ndigits)
{
  //set fixed precision
  ostringstream out;
  out<<std::fixed;
  
  //get the negative of pow of 10 needed to make all syst have at least 1 digit
  double mel=*max_element(errors.begin(),errors.end());
  if(mel==0) out<<ave<<"(0)";
  else
    if(isnan(mel)) out<<ave<<"(nan)";
    else
      {
	int lem=(int)ceil(log(mel)/M_LN10);
	int digit_to_round=lem-ndigits;
	out.precision((lem<ndigits)*max(-digit_to_round,0));
	
	//truncate
	out<<round_to_digit(ave,digit_to_round);
	if(digit_to_round<=-ndigits) out.precision(0);
	for(auto err : errors) out<<"("<<shift_by_10(round_to_digit(err,digit_to_round),(lem<=0)*(ndigits-lem))<<")";
      }
  
  return out.str();
}
