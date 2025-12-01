#include "timer.hpp"
#include <map>
#include <string>
#include <memory>
extern "C" {

  typedef void* timer_handle;

  timer_handle mytimer_create(const char* name) {
      std::string sname(name);
      return new CTimer(sname);
  }

  void mytimer_destroy(timer_handle handle) {
      delete static_cast<CTimer*>(handle);
  }

  void mytimer_print() {
      CTimer::print_timing_results();
  }

}
