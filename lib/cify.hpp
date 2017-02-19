#ifndef _CIFY_HPP
#define _CIFY_HPP

#include <utility>

namespace
{
template <typename F, int I, typename L, typename R, typename ...A>
inline F cify(L&& l, R (*)(A...) noexcept(noexcept(
  ::std::declval<F>()(::std::declval<A>()...))))
  {
    static thread_local L l_(::std::forward<L>(l));
    static thread_local bool full;
    
    if(full)
      {
	l_.~L();
	
	new (static_cast<void*>(&l_)) L(::std::forward<L>(l));
      }
    else full=true;
    
    struct S
    {
      static R f(A... args) noexcept(noexcept(::std::declval<F>()(::std::forward<A>(args)...)))
      {return l_(::std::forward<A>(args)...);}
    };
    
    return &S::f;
  }
}

template <typename F,int I=0,typename L>
F cify(L&& l)
{return cify<F, I>(::std::forward<L>(l),F());}


#endif
