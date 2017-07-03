#ifndef _STOPWATCH_HPP
#define _STOPWATCH_HPP

#include <chrono>
#include <string>

#include <tools.hpp>

using namespace std;

//! measure time
using instant_t=chrono::time_point<chrono::steady_clock>;
inline instant_t take_time()
{return chrono::steady_clock::now();}

//! difference between two measures
using duration_t=decltype(take_time()-take_time());

//! implment a stopwatch
class stopwatch_t : duration_t
{
  //! name of the stopwatch
  const string name;
  
  //! number of measurements
  size_t nmeas;
  
  //! last measurement
  instant_t last;
  
  //! true if the stopwatch has been started
  bool started;
  
public:
  stopwatch_t(const string &name) : name(name) {reset();}
  
  //! descritpion
  string descr() const
  {return name;}
  
  //! reference to base type
  duration_t& base()
  {return (duration_t&)(*this);}
  
  //! const reference to base type
  const duration_t& base() const
  {return (const duration_t&)(*this);}
  
  //! start the stopwatch
  void start()
  {
    if(started==true) CRASH("Trying to start an already started stopwatch!");
    last=take_time();
    started=true;
  }
  
  //! stop the stopwatch
  void stop()
  {
    if(started==false) CRASH("Trying to stop a non-started stopwatch!");
    base()+=take_time()-last;
    nmeas++;
    started=false;
  }
  
  //! reset the timewatch
  void reset()
  {
    nmeas=0;
    started=false;
    last=take_time();
    base()=last-last;
  }
  
  //! return the total time
  duration_t tot() const
  {
    if(started) CRASH("Tryng to ask for the total time of a running stopwatch!");
    return base();
  }
  
  //! return the average
  duration_t ave() const
  {return tot()/nmeas;}
  
  //! return the number of measures
  size_t count() const
  {return nmeas;}
  
  //! return the sum of two stopwatch
  stopwatch_t operator+=(const stopwatch_t &oth)
  {
    if(name!=oth.name) CRASH("Trying to combine two stopwatch with different name, %s and %s",name.c_str(),oth.name.c_str());
    
    nmeas+=oth.nmeas;
    base()+=oth.base();
    
    return *this;
  }
};

//! print elapsed time
inline ostream& operator<<(ostream &out,const duration_t &diff)
{
  double el_nano=chrono::duration<double,nano>(diff).count();
  if(el_nano<1000) return out<<el_nano<<" ns";
  
  double el_micro=chrono::duration<double,micro>(diff).count();
  if(el_micro<1000) return out<<el_micro<<" us";
  
  double el_milli=chrono::duration<double,milli>(diff).count();
  if(el_milli<1000) return out<<el_milli<<" ms";
  
  double el_sec=chrono::duration<double>(diff).count();
  return out<<el_sec<<" s";
}

//! return elapsed time from a given moment
inline duration_t elapsed_time(const instant_t &start)
{return take_time()-start;}

//! print a stopwatch
inline ostream& operator<<(ostream &out,const stopwatch_t &s)
{return out<<"Time to "<<s.descr()<<" "<<s.count()<<" time: "<<s.tot();}

//! ratio of two stopwatches
class stopwatch_ratio_t
{
public:
  //! numerator
  const stopwatch_t &num;
  
  //! denominator
  const stopwatch_t &den;
  
  stopwatch_ratio_t(const stopwatch_t &num,const stopwatch_t &den) : num(num),den(den) {}
};

//! print a stopwatch
inline ostream& operator<<(ostream &out,const stopwatch_ratio_t &s)
{return out<<s.num<<",\t relative to "<<s.den.descr()<<":\t "<<s.num.tot()*100.0/s.den.tot()<<" %";}

//! return a stopwatch_ratio to write the used time in proportion of another one
inline stopwatch_ratio_t operator/(const stopwatch_t &num,const stopwatch_t &den)
{return stopwatch_ratio_t(num,den);}

//! reduction for stopwatch
#pragma omp declare reduction(+: stopwatch_t : omp_out+=omp_in) initializer(omp_priv=stopwatch_t(omp_orig.descr()))

#endif
