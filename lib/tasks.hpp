#ifndef _TASKS_HPP
#define _TASKS_HPP

#include <future>
#include <memory>
#include <vector>

using namespace std;

using incapsulated_task_t=function<void(bool)>;

const bool RECYCLE=true,DONT_RECYCLE=false;

//! single task
template<typename F,typename... Rest>
incapsulated_task_t* incapsulate_task(F &&f,Rest&&... rest)
{
  auto pck=make_shared<packaged_task<decltype(f(rest...))()>>(bind(forward<F>(f),forward<Rest>(rest)...));
  return new function<void(bool)>([pck](bool recycle){(*pck)();if(recycle) pck->reset();});
}

//! list of tasks
class task_list_t : public vector<incapsulated_task_t*>
{
public:
  void assolve_all(bool recycle=DONT_RECYCLE)
  {
#pragma omp parallel for
  for(size_t itask=0;itask<this->size();itask++)
    (*((*this)[itask]))(recycle);
  }
};

#endif
