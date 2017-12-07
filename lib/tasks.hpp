#ifndef _TASKS_HPP
#define _TASKS_HPP

#include <future>
#include <memory>
#include <vector>

using namespace std;

//! incapsulate a task to be used at later stage
using incapsulated_task_t=function<void(bool)>;

//! allows the task to be called again
const bool RECYCLE=true;

//! don't allow the task to be called again
const bool DONT_RECYCLE=false;

//! single task
template<typename F,typename... Rest>
incapsulated_task_t incapsulate_task(F &&f,Rest&&... rest)
{
  auto pck=make_shared<packaged_task<decltype(f(rest...))()>>
    (bind(forward<F>(f),forward<Rest>(rest)...));
  return function<void(bool)>
    ([pck](bool recycle)
     {
       (*pck)();
       if(recycle) pck->reset();
     });
}

//! list of tasks
class task_list_t : public vector<incapsulated_task_t>
{
public:
  void assolve_all(const bool recycle=DONT_RECYCLE)
  {
#pragma omp parallel for
  for(size_t itask=0;itask<this->size();itask++)
    ((*this)[itask])(recycle);
  }
};

#endif
